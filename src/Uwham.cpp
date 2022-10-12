#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const WhamInput& input)
:Wham(input)
{
    whamPack_->Readbool("BAR", ParameterPack::KeyType::Optional, BAR_);
    bool readE=whamPack_->Readbool("ErrorAnalysis", ParameterPack::KeyType::Optional, Error_);
    if (readE){
        whamPack_->ReadNumber("ErrorIteration", ParameterPack::KeyType::Required, ErrorIter_);
    }

    registerOutput("normalization", [this](std::string name)-> void {this->printNormalization(name);});
    registerOutput("pji", [this](std::string name)->void{this->printPji(name);});
    registerOutput("lnwji", [this](std::string name)->void{this->printlnwji(name);});
    registerOutput("derivative", [this](std::string name) -> void{this -> printderivative(name);});
    registerOutput("reweightFE", [this](std::string name) -> void {this -> printReweightFE(name);});
    registerOutput("KL_divergence", [this](std::string name) -> void {this -> printKL(name);});
    registerOutput("FE_dim", [this](std::string name) -> void {this -> printFEdim(name);});

    // check if the outputs are registered
    isRegistered();

    // construct the BUki matrix 
    initializeBUki();

    // make BAR initial guess
    MakeInitialGuess(BUki_, N_, fk_);
 
    // Now read in what type of calculation are you letting Uwham do
    initializeStrat(BUki_, N_, strategies_);

    // bin the data upfront
    bindata(xi_, MapBinIndexTolnwjiIndex_, DataBinIndex_);
}


void Uwham::initializeBUki()
{
    std::cout << "The total data in BUki is " << Biases_.size() * xi_.size() << std::endl;
    std::cout << "The space it takes is around " << Biases_.size() * xi_.size() * sizeof(Real) / 1e6 << " Mb" << std::endl;
    calculateBUki(xi_, BUki_); 
}

void Uwham::calculateBUki(const std::vector<std::vector<Real>>& xi, Matrix<Real>& BUki)
{
    BUki.resize(Biases_.size(), xi.size());
    #pragma omp parallel for
    for (int i=0;i<xi.size();i++)
    {
        for (int j=0;j<Biases_.size();j++){
            Real val = Biases_[j]->calculate(xi[i]); 
            BUki(j,i) = Biases_[j]->getBeta()*val;
        }
    }
}

void Uwham::MakeInitialGuess(const Matrix<Real>& BUki, const std::vector<Real>& N, std::vector<Real>& fk)
{
    // inital guess for fk 
    fk.clear();
    fk.resize(N.size(),0.0);

    // check if we are doing Bennet Acceptance Ratio (BAR) for initial guess 
    if (BAR_){
        // map from group index to point index in xi
        std::vector<std::vector<int>> GroupIndex;
        MakeGroupPointMap(N, GroupIndex);

        fk = std::vector<Real>(N.size(),0.0);

        for (int i=0;i<GroupIndex.size()-1;i++)
        {
            int forwardSize = GroupIndex[i].size();
            int backwardSize= GroupIndex[i+1].size();

            // forward work
            std::vector<Real> w_F(forwardSize);
            std::vector<Real> w_B(backwardSize);

            int k = i;
            int l = i+1;

            for (int j=0;j<forwardSize;j++){
                w_F[j] = BUki(l,GroupIndex[k][j]) - BUki(k, GroupIndex[k][j]);
            }

            for (int j=0;j<backwardSize;j++){
                w_B[j] = BUki(k, GroupIndex[l][j]) - BUki(l, GroupIndex[l][j]);
            }

            Real DeltaF = WhamTools::CalculateDeltaFBarIterative(w_F, w_B);

            fk[l] = fk[k] + DeltaF;
        }
    }
    else
    {
        fk = std::vector<Real>(N.size(),0.0);
    }
}

void Uwham::MakeGroupPointMap(const std::vector<Real>& N, std::vector<std::vector<int>>& GroupIndex)
{
    GroupIndex.clear();
    GroupIndex.resize(N_.size());

    int initial=0;
    for (int i=0;i<N.size();i++)
    {
        std::vector<int> temp(N[i],0);
        std::iota(temp.begin(), temp.end(), initial);
        GroupIndex[i] = temp;
        initial = initial + N[i];
    }
}

void Uwham::initializeStrat(Matrix<Real>& BUki, std::vector<Real>& N, std::vector<stratptr>& strategies)
{
    strategies.clear();

    auto whampack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    std::map<std::string, UWhamCalculationStrategy*> MapNameToStrat;
    std::vector<std::string> strategyNames;

    auto stratPacks = whampack->findParamPacks("Uwhamstrategy", ParameterPack::KeyType::Required);

    for (auto s : stratPacks)
    {
        std::string strattype;
        std::string name;
        UwhamStrategyInput input = {BUki, N, const_cast<ParameterPack&>(*s)};

        s -> ReadString("type", ParameterPack::KeyType::Required, strattype);
        auto sptr = UwhamCalculationStrategyRegistry::Factory::instance().create(strattype, input);
        MapNameToStrat.insert(std::make_pair(sptr -> getName(), sptr));
    }

    // read a vector of string that represents the order of optimization that we want to do , usually LBFGS --> adaptive     
    whampack->ReadVectorString("strategyNames", ParameterPack::KeyType::Required, strategyNames);

    for (auto s : strategyNames)
    {
        auto stratit = MapNameToStrat.find(s);
        ASSERT((stratit != MapNameToStrat.end()), "Strategy name " << s << " not found.");

        // transfer ownership of the pointer
        strategies.push_back(stratptr(stratit->second));
    }
}

void Uwham::bindata(std::vector<std::vector<Real>>& xi, std::map<std::vector<int>, std::vector<int>>& map, std::vector<std::vector<int>>& DataBinIndex)
{
    // clear the bin index for each of the data point
    DataBinIndex.clear();
    DataBinIndex.resize(xi.size());

    // clear the map
    map.clear();

    OpenMP::OpenMP_buffer<std::map<std::vector<int>, std::vector<int>>> mapBuffer;
    mapBuffer.set_master_object(map);

    #pragma omp parallel
    {
        // access the local lnwji index map
        auto& localmap  = mapBuffer.access_buffer_by_id();

        // clear the local maps
        localmap.clear();

        #pragma omp for
        for (int i=0;i<xi.size();i++){
            std::vector<Real>& x = xi[i];

            std::vector<int> BinIndex;
            BinIndex.resize(Bins_.size());

            bool isInRange = true; 

            for (int j=0;j<Bins_.size();j++){
                // User input is 1 based
                int dim = Bins_[j].getDimension() - 1;

                // If this data point is out of range for one of the bins, it is out of range, so we break
                if (! Bins_[j].isInRange(xi[i][dim])){
                    // break from the for loop
                    isInRange = false;
                    break;
                }
            
                int index = Bins_[j].findBin(xi[i][dim]);
                BinIndex[j] = index;
            }

            if (isInRange){
                // add the indices to the binned data vector
                templatetools::InsertIntoVectorMap(BinIndex, i, localmap);
                DataBinIndex[i] = BinIndex;
            }
        }
    }

    // Perform the same operation but for map bin index to ln index
    for (auto m = mapBuffer.beginworker(); m != mapBuffer.endworker(); m ++)
    {
        for (auto it = m->begin(); it != m -> end(); it++)
        {
            auto Indexit = map.find(it -> first);

            // if the original map also has the bin index
            if (Indexit != map.end())
            {
                Indexit -> second.insert(Indexit -> second.end(), it -> second.begin(), it ->second.end());
            }
            else
            {
                map.insert(std::make_pair(it -> first, it -> second));
            }
        }
    }
}

void Uwham::calculateFreeEnergy(const std::vector<Real>& lnwji, std::map<std::vector<int>, std::vector<int>>& map, std::map<std::vector<int>, Real>& FE)
{
    FE.clear();

    for (auto it = map.begin();it != map.end();it++)
    {
        auto& l = it -> second;
        std::vector<Real> ones(l.size(), 1.0);
        std::vector<Real> lnwji_bin(l.size(),0.0);

        for (int i=0;i<l.size();i++)
        {
            lnwji_bin[i] = lnwji[l[i]];
        }

        Real wji = WhamTools::LogSumExp(lnwji_bin, ones);

        FE.insert(std::make_pair(it->first, wji));
    }

    return;
}

void Uwham::calculate()
{
    // Start calculation --> first calculate using the entire data set 
    for (int i=0;i<strategies_.size();i++){
        strategies_[i] -> calculate(fk_);
        fk_ = strategies_[i] -> getFk_();
        lnwji_ = strategies_[i] -> getlnwji_();
    }

    // calculate the free energy 
    calculateFreeEnergy(lnwji_, MapBinIndexTolnwjiIndex_, FreeEnergy_);

    // We should also reweight the data to each of the simulations
    // resize lnpji to the size of 'Number of biases'
    lnpji_.resize(BUki_.getNR(), std::vector<Real>(lnwji_.size(),0.0));
    reweightFE_.resize(BUki_.getNR());

    // Reweight to various biases 
    for (int i=0;i<BUki_.getNR();i++)
    {
        #pragma omp parallel for
        for (int j=0;j<lnwji_.size();j++)
        {
            lnpji_[i][j] = fk_[i] - BUki_(i,j) + lnwji_[j];
        }
    }

    // only do this procedure if we are not doing combined input
    if (! combined_input_)
    {
        for (int i=0;i<BUki_.getNR();i++)
        {
            for (auto it = MapBinIndexTolnwjiIndex_.begin(); it != MapBinIndexTolnwjiIndex_.end(); it ++)
            {
                auto& l  = it -> second;
                std::vector<Real> lnpi;

                for (int j=0;j< l.size(); j++)
                {
                    lnpi.push_back(lnpji_[i][l[j]]);
                }

                std::vector<Real> ones(l.size(), 1.0);
                Real pi = -WhamTools::LogSumExp(lnpi, ones);

                reweightFE_[i].insert(std::make_pair(it -> first, pi));
            }
        }

        // Now let's calculate KL divergence
        KL_divergence_.resize(BUki_.getNR(),0.0);
        for (int i=0;i<BUki_.getNR();i++)
        {
            for (auto it = dataFE_[i].begin(); it != dataFE_[i].end(); it ++)
            {
                auto Index = it -> first;

                Real ref_val = it -> second;
                Real prob = std::exp(-ref_val);
                Real val = reweightFE_[i].find(Index) -> second;

                KL_divergence_[i] = prob * (-ref_val + val);
            }
        }
    }

    // get the FE in each of the dimensions
    ReduceFEDimension();

    // if we want to calculate error 
    if (Error_){
        calculateError();
    }
}

void Uwham::calculateError()
{
    for (int i=0;i<ErrorIter_;i++)
    {
        std::vector<std::vector<Real>> X;
        std::vector<Real> N;

        for (int i=0;i<VectorTimeSeries_.size();i++)
        {
            auto& ts = VectorTimeSeries_[i];
            std::vector<std::vector<Real>> sample = ts->getIndependentsample();
            X.insert(X.end(), sample.begin(), sample.end());
            N.push_back(sample.size());
        }

        // calculate BUki
        Matrix<Real> BUki;
        calculateBUki(X, BUki);

        // bin the data 
        std::map<std::vector<int>,std::vector<int>> map;
        std::vector<std::vector<int>> dbinIndex;
        bindata(X, map, dbinIndex);

        // make initial guess for fk = -log(Qk)
        std::vector<Real> fk_guess;
        std::vector<Real> lnwji;
        MakeInitialGuess(BUki, N, fk_guess);

        // make new strategies
        std::vector<stratptr> strategies;
        std::map<std::vector<int>, Real> FE;
        initializeStrat(BUki, N, strategies);

        for (int i=0;i<strategies.size();i++)
        {
            strategies[i] -> calculate(fk_guess);
            fk_guess = strategies[i]->getFk_();
            lnwji = strategies[i]->getlnwji_();
        }

        // calculate the Free Energy
        calculateFreeEnergy(lnwji, map, FE);

        for (auto it = FE.begin(); it != FE.end(); it ++)
        {
            std::vector<int> copyK = it -> first;
            templatetools::InsertIntoVectorMap(copyK, it -> second, ErrorFEMap_);
        }
    }

    for (auto it = ErrorFEMap_.begin(); it != ErrorFEMap_.end(); it ++)
    {
        Real mean=0.0;
        Real var =0.0;
        for (int i=0;i<it->second.size();i++)
        {
            mean += it -> second[i];
        }

        mean = mean / ErrorIter_;
        for (int i=0;i<it->second.size();i++)
        {
            var += std::pow((it->second[i] - mean),2.0);
        }

        var = var + (ErrorIter_ - it->second.size()) * (mean*mean);
        var = var / ErrorIter_;

        Real std = std::sqrt(var);

        ErrorMap_.insert(std::make_pair(it -> first, std));
        MeanMap_.insert(std::make_pair(it->first, mean));
    }
}

void Uwham::ReduceFEDimension()
{
    FE_dim_.clear();
    FE_dim_.resize(Bins_.size());

    for (int i=0;i<Bins_.size();i++)
    {
        // map each of the bin num in various dimension to lnwji index 
        std::map<int, std::vector<int>> MapDimBinNumTolnwjiIndex;
        for (auto it = MapBinIndexTolnwjiIndex_.begin(); it != MapBinIndexTolnwjiIndex_.end();it++)
        {
            int binnum = it -> first[i];
            auto value = it -> second;

            auto itt = MapDimBinNumTolnwjiIndex.find(binnum);

            if (itt != MapDimBinNumTolnwjiIndex.end())
            {
                itt -> second.insert(itt -> second.end(), value.begin(), value.end());
            }
            else
            {
                MapDimBinNumTolnwjiIndex.insert(std::make_pair(binnum, value));
            }
        }
        // Now using those index, find the -log(sum(exp(lnwji)))
        for (auto it = MapDimBinNumTolnwjiIndex.begin(); it != MapDimBinNumTolnwjiIndex.end();it++)
        {
            std::vector<Real> lnwji_dim;

            for (int j=0;j<it -> second.size();j++)
            {
                lnwji_dim.push_back(lnwji_[it ->second[j]]);
            }

            std::vector<Real> ones(lnwji_dim.size(),1.0);
            Real Fe = -1.0 * WhamTools::LogSumExpOMP(lnwji_dim, ones);

            FE_dim_[i].insert(std::make_pair(it -> first, Fe));
        }
    }
}

int Uwham::getNumBinsPerDimension(int num)
{
    ASSERT((num < Bins_.size()), "The dimension provided is larger than the total number of dimensions which is " << Bins_.size());

    return Bins_[num].getNumbins();
}

                            /*************************************
                             ********** print functions **********
                             ************************************/  
void Uwham::printOutput()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }
}

void Uwham::printlnwji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<lnwji_.size();i++)
    {
        ofs << lnwji_[i] << std::endl;
    }

    ofs.close();
}

void Uwham::printPji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);
    ofs << "#";

    for (int i=0;i<Bins_.size();i++)
    {
        ofs << "Position" << Bins_[i].getDimension() << "\t";
    }

    for (int i=0;i<Bins_.size();i++)
    {
        ofs << "Index" << Bins_[i].getDimension() << "\t";
    }

    ofs << "pji\tF";

    if (Error_)
    {
        ofs << "\tStd";
        ofs << "\tMean";
    }

    ofs << "\n";

    for (auto it = FreeEnergy_.begin();it != FreeEnergy_.end();it++)
    {
        auto& index = it -> first;
        for (int i=0;i<Bins_.size();i++)
        {
            Real pos = Bins_[i].getLocationOfBin(index[i]);
            ofs << pos << " ";
        }

        for (int i=0;i<Bins_.size();i++)
        {
            ofs << index[i] << " ";
        }

        ofs << it -> second << " ";
        ofs << (-1.0)*(it -> second); 

        if (Error_)
        {
            auto it = ErrorMap_.find(index);

            if (it != ErrorMap_.end())
            {
                ofs << " " << it -> second;
            }
            else
            {
                ofs << " " << 0;
            }

            auto itM = MeanMap_.find(index);

            if ( itM != MeanMap_.end())
            {
                ofs << " " << itM -> second;
            }
            else
            {
                ofs << " " << 0;
            }
        }

        ofs << "\n";
    }
    ofs.close();
}

void Uwham::printNormalization(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int Nsim = BUki_.getNR();
    ofs << std::fixed << std::setprecision(precision_);
    ofs << "# normalization constants" << "\n";

    for (int i=0;i<Nsim;i++)
    {
        ofs << fk_[i] << "\n";
    }
    ofs.close();
}

void Uwham::printKL(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<KL_divergence_.size();i++)
    {
        ofs << KL_divergence_[i] << "\n";
    }
    ofs.close();
}

void Uwham::printReweightFE(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<reweightFE_.size();i++)
    {
        for (auto it = reweightFE_[i].begin(); it != reweightFE_[i].end(); it ++)
        {
            ofs << i+1 << " ";
            for (int j=0;j<it->first.size();j++)
            {
                ofs << it -> first[j] << " ";
            }
            ofs << it -> second << "\n";
        }
    }

    ofs.close();
}

void Uwham::printFEdim(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    ofs << "# dimension bin FE \n";
    for (int i=0;i<FE_dim_.size();i++)
    {
        for (auto it = FE_dim_[i].begin(); it != FE_dim_[i].end(); it ++)
        {
            ofs << i+1 << " " << it -> first << " " << it -> second << "\n";
        }
    }

    ofs.close();
}

void Uwham::printderivative(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    // we have to calculate the derivative
    int dim = Bins_.size();

    std::vector<std::vector<Real>> dFdN(Biases_.size(), std::vector<Real>(dim,0.0));
    std::vector<std::vector<Real>> avgs(Biases_.size(), std::vector<Real>(dim,0.0));

    // each bias has a point
    for (int i=0;i<Biases_.size();i++)
    {
        std::vector<Real> avg(dim,0.0);
        for (int j=0;j<xi_.size();j++)
        {
            Real factor = fk_[i] - BUki_(i,j) + lnwji_[j];
            factor = std::exp(factor);

            for (int k=0;k<dim;k++)
            {
                avg[k] += factor * xi_[j][k];
            }
        }

        avgs[i] = avg;
        std::vector<Real> dF = Biases_[i] -> calculateForce(avg);
        ASSERT((dF.size() == dim), "The dimensions don't match.");

        dFdN[i] = dF;
    }

    for (int i=0;i<Biases_.size();i++)
    {
        for (int j=0;j<dim;j++)
        {
            ofs << avgs[i][j] << " ";
        }

        for(int j=0;j<dim;j++)
        {
            ofs << dFdN[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}