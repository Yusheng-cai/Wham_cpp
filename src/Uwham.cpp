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

    // check if the outputs are registered
    isRegistered();

    // construct the BUki matrix 
    initializeBUki();

    // make BAR initial guess
    BARInitialGuess();
 
    // Now read in what type of calculation are you letting Uwham do
    initializeStrat();
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
    BUki_.resize(Biases_.size(), xi_.size());
    #pragma omp parallel for
    for (int i=0;i<xi_.size();i++)
    {
        for (int j=0;j<Biases_.size();j++)
        {
            Real val = Biases_[j]->calculate(xi_[i]); 
            BUki_(j,i) = Biases_[j]->getBeta()*val;
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
    std::cout << "Calculation took " << diff.count() << std::endl;

    binneddata_.resize(xi_.size());

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
            auto it = MapBinIndexToVectorlnwji_.find(BinIndex);
            auto it2= MapBinIndexTolnwjiIndex_.find(BinIndex);

            if ( it == MapBinIndexToVectorlnwji_.end())
            {
                std::vector<Real> lnwji_vec_;
                std::vector<int> lnwjiIndex_vec_;

                lnwji_vec_.push_back(lnwji_[i]);
                lnwjiIndex_vec_.push_back(i);

                MapBinIndexToVectorlnwji_.insert(std::make_pair(BinIndex, lnwji_vec_));
                MapBinIndexTolnwjiIndex_.insert(std::make_pair(BinIndex, lnwjiIndex_vec_));
            }
            else
            {
                it -> second.push_back(lnwji_[i]);
                it2 -> second.push_back(i);
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