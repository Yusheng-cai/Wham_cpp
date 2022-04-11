#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const WhamInput& input)
:Wham(input)
{
    bool real = whamPack_->Readbool("BAR", ParameterPack::KeyType::Optional, BAR_);

    registerOutput("normalization", [this](std::string name)-> void {this->printNormalization(name);});
    registerOutput("pji", [this](std::string name)->void{this->printPji(name);});
    registerOutput("lnwji", [this](std::string name)->void{this->printlnwji(name);});
    registerOutput("derivative", [this](std::string name) -> void{this -> printderivative(name);});
    registerOutput("derivativeNormTS", [this](std::string name) -> void {this -> printderivativeNormTS(name);});
    registerOutput("reweightFE", [this](std::string name) -> void {this -> printReweightFE(name);});
    registerOutput("KL_divergence", [this](std::string name) -> void {this -> printKL(name);});
    registerOutput("FE_dim", [this](std::string name) -> void {this -> printFEdim(name);});

    // check if the outputs are registered
    isRegistered();

    // construct the BUki matrix 
    initializeBUki();

    // make BAR initial guess
    BARInitialGuess();
 
    // Now read in what type of calculation are you letting Uwham do
    initializeStrat();
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

void Uwham::printderivativeNormTS(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    const auto& derivativeNorm = strat_->getNorm();

    ofs << " Iteration \t DerivativeNorm \n";
    for (int i=0;i<derivativeNorm.size();i++)
    {
        ofs << i+1 << "\t" << derivativeNorm[i] << "\n";
    }
    ofs.close();
}

void Uwham::printderivative(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    // we have to calculate the derivative
    const auto& lnwji = strat_ -> getlnwji_();
    const auto& fk    = strat_ -> getFk_();
    int dim = Bins_.size();

    std::vector<std::vector<Real>> dFdN(Biases_.size(), std::vector<Real>(dim,0.0));
    std::vector<std::vector<Real>> avgs(Biases_.size(), std::vector<Real>(dim,0.0));

    // each bias has a point
    for (int i=0;i<Biases_.size();i++)
    {
        std::vector<Real> avg(dim,0.0);
        for (int j=0;j<xi_.size();j++)
        {
            Real factor = fk[i] - BUki_(i,j) + lnwji[j];
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
        for (int j=0;j<Biases_.size();j++)
        {
            Real val = Biases_[j]->calculate(xi[i]); 
            BUki(j,i) = Biases_[j]->getBeta()*val;
        }
    }
}

void Uwham::BARInitialGuess()
{
    if (BAR_)
    {
        // map from group index to point index in xi
        MakeGroupPointMap();

        fk_ = std::vector<Real>(N_.size(),0.0);

        for (int i=0;i<GroupIndex_.size()-1;i++)
        {
            int forwardSize = GroupIndex_[i].size();
            int backwardSize= GroupIndex_[i+1].size();

            // forward work
            std::vector<Real> w_F(forwardSize);
            std::vector<Real> w_B(backwardSize);

            int k = i;
            int l = i+1;

            for (int j=0;j<forwardSize;j++)
            {
                w_F[j] = BUki_(l,GroupIndex_[k][j]) - BUki_(k, GroupIndex_[k][j]);
            }

            for (int j=0;j<backwardSize;j++)
            {
                w_B[j] = BUki_(k, GroupIndex_[l][j]) - BUki_(l, GroupIndex_[l][j]);
            }

            Real DeltaF = WhamTools::CalculateDeltaFBarIterative(w_F, w_B);

            fk_[l] = fk_[k] + DeltaF;
            std::cout << "fk " << l << " = " << fk_[l] << "\n";
        }
    }
    else
    {
        std::cout << "initializing fk to zeros." << "\n";
        fk_ = std::vector<Real>(N_.size(),0.0);
    }

    whamPack_->ReadVectorNumber("fk", ParameterPack::KeyType::Optional, fk_);
}

void Uwham::MakeGroupPointMap()
{
    GroupIndex_.clear();
    GroupIndex_.resize(N_.size());

    int initial=0;
    for (int i=0;i<N_.size();i++)
    {
        std::vector<int> temp(N_[i],0);
        std::iota(temp.begin(), temp.end(), initial);
        GroupIndex_[i] = temp;
        initial = initial + N_[i];
    }
}


void Uwham::initializeStrat()
{
    auto whampack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    auto stratPacks = whampack->findParamPacks("Uwhamstrategy", ParameterPack::KeyType::Required);
    for (auto s : stratPacks)
    {
        std::string strattype;
        std::string name;
        UwhamStrategyInput input = {BUki_, N_, const_cast<ParameterPack&>(*s), fk_};
        int index = strategies_.size();
        s -> ReadString("type", ParameterPack::KeyType::Required, strattype);
        strategies_.push_back(stratptr(UwhamCalculationStrategyRegistry::Factory::instance().create(strattype, input)));
        MapNameToStrat_.insert(std::make_pair(strategies_[index]->getName(), strategies_[index].get()));
    }


    // read a vector of string that represents the order of optimization that we want to do , usually LBFGS --> adaptive     
    whampack->ReadVectorString("strategyNames", ParameterPack::KeyType::Required, strategyNames_);

    for (auto s : strategyNames_)
    {
        auto stratit = MapNameToStrat_.find(s);
        ASSERT((stratit != MapNameToStrat_.end()), "Strategy name " << s << " not found.");
    }
}

void Uwham::calculate()
{
    // echo number of threads we are using 
    int numthreads = 0;
    #pragma omp parallel
    {
        #pragma omp critical
        numthreads += 1;
    }

    // Start calculation
    std::cout << "We are using " << numthreads << " OpenMP threads for this operation." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto s : strategyNames_)
    {
        auto strat = MapNameToStrat_.find(s)->second;
        strat -> setFk(fk_);
        strat -> calculate();
        fk_ = strat -> getFk_();
        lnwji_ = strat -> getlnwji_();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Calculation took " << diff.count() << " us." << std::endl;

    // Now we perform the binning 
    binneddata_.resize(xi_.size());
    auto s = std::chrono::high_resolution_clock::now();
    MapBinIndexToVectorlnwjiBuffer_.set_master_object(MapBinIndexToVectorlnwji_);
    MapBinIndexTolnwjiIndexBuffer_.set_master_object(MapBinIndexTolnwjiIndex_);
    #pragma omp parallel
    {
        // access the local vector lnwji map 
        auto& localVectorlnwjiMap = MapBinIndexToVectorlnwjiBuffer_.access_buffer_by_id();
        // access the local lnwji index map
        auto& locallnwjiIndexMap  = MapBinIndexTolnwjiIndexBuffer_.access_buffer_by_id();

        // clear the local maps
        localVectorlnwjiMap.clear();
        locallnwjiIndexMap.clear();

        #pragma omp for
        for (int i=0;i<xi_.size();i++)
        {
            auto& x = xi_[i];
            ASSERT((x.size() == Bins_.size()), "The dimension of the data = " << x.size() << " while the dimension of bins = " << Bins_.size());
            std::vector<int> BinIndex;
            BinIndex.resize(Bins_.size());

            bool isInRange = true; 

            for (int j=0;j<Bins_.size();j++)
            {
                // User input is 1 based
                int dim = Bins_[j].getDimension() - 1;

                // If this data point is out of range for one of the bins, it is out of range, so we break
                if (! Bins_[j].isInRange(xi_[i][dim]))
                {
                    // break from the for loop
                    isInRange = false;
                    break;
                }
            
                int index = Bins_[j].findBin(xi_[i][dim]);
                BinIndex[j] = index;
            }

            if (isInRange)
            {
                // add the indices to the binned data vector
                binneddata_[i] = BinIndex;
                auto it = localVectorlnwjiMap.find(BinIndex);
                auto it2= locallnwjiIndexMap.find(BinIndex);

                if ( it == localVectorlnwjiMap.end())
                {
                    std::vector<Real> lnwji_vec_;
                    std::vector<int> lnwjiIndex_vec_;

                    lnwji_vec_.push_back(lnwji_[i]);
                    lnwjiIndex_vec_.push_back(i);

                    localVectorlnwjiMap.insert(std::make_pair(BinIndex, lnwji_vec_));
                    locallnwjiIndexMap.insert(std::make_pair(BinIndex, lnwjiIndex_vec_));
                }
                else
                {
                    it -> second.push_back(lnwji_[i]);
                    it2 -> second.push_back(i);
                }
            }
        }
    }

    // Now we combine all the maps for map bin index to vector lnwji
    for (auto map = MapBinIndexToVectorlnwjiBuffer_.beginworker(); map != MapBinIndexToVectorlnwjiBuffer_.endworker(); map ++)
    {
        for (auto it = map->begin(); it != map -> end(); it++)
        {
            auto Vectorit = MapBinIndexToVectorlnwji_.find(it ->first);

            // if the original map also has the bin index
            if (Vectorit != MapBinIndexToVectorlnwji_.end())
            {
                Vectorit -> second.insert(Vectorit -> second.end(), it -> second.begin(), it ->second.end());
            }
            else
            {
                MapBinIndexToVectorlnwji_.insert(std::make_pair(it -> first, it -> second));
            }
        }
    }

    // Perform the same operation but for map bin index to ln index
    for (auto map = MapBinIndexTolnwjiIndexBuffer_.beginworker(); map != MapBinIndexTolnwjiIndexBuffer_.endworker(); map ++)
    {
        for (auto it = map->begin(); it != map -> end(); it++)
        {
            auto Indexit = MapBinIndexTolnwjiIndex_.find(it -> first);

            // if the original map also has the bin index
            if (Indexit != MapBinIndexTolnwjiIndex_.end())
            {
                Indexit -> second.insert(Indexit -> second.end(), it -> second.begin(), it ->second.end());
            }
            else
            {
                MapBinIndexTolnwjiIndex_.insert(std::make_pair(it -> first, it -> second));
            }
        }
    }

    for (auto it = MapBinIndexToVectorlnwji_.begin();it != MapBinIndexToVectorlnwji_.end();it++)
    {
        auto& l = it -> second;
        std::vector<Real> ones(l.size(), 1.0);

        Real wji = WhamTools::LogSumExp(l, ones);

        MapBinIndexToWji_.insert(std::make_pair(it->first, wji));
    }

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

    // get the FE in each of the dimensions
    ReduceFEDimension();

    auto e = std::chrono::high_resolution_clock::now();
    auto d = std::chrono::duration_cast<std::chrono::microseconds>(e - s);
    std::cout << "Binning took " << d.count() << " us." << "\n";
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

int Uwham::getNumBinsPerDimension(int num)
{
    ASSERT((num < Bins_.size()), "The dimension provided is larger than the total number of dimensions which is " << Bins_.size());

    return Bins_[num].getNumbins();
}

void Uwham::printlnwji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    const auto& lnwji = strat_ -> getlnwji_();

    for (int i=0;i<lnwji.size();i++)
    {
        ofs << lnwji[i] << std::endl;
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
    ofs << "\n";

    for (auto it = MapBinIndexToWji_.begin();it != MapBinIndexToWji_.end();it++)
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
        ofs << (-1.0)*(it -> second) << "\n"; 
    }
    ofs.close();
}

void Uwham::printOutput()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }
}

void Uwham::finishCalculate()
{
}