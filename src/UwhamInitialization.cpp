#include "Uwham.h"

#include <numeric>

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
    registerOutput("ErrorFE", [this](std::string name) -> void {this->printErrroFE(name);});

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
    if (verbose_){
        std::cout << "The total data in BUki is " << Biases_.size() * xi_.size() << std::endl;
        std::cout << "The space it takes is around " << Biases_.size() * xi_.size() * sizeof(Real) / 1e6 << " Mb" << std::endl;
    }
    calculateBUki(xi_, BUki_);
}

void Uwham::calculateBUki(const std::vector<std::vector<Real>>& xi, Matrix<Real>& BUki)
{
    BUki.resize(Biases_.size(), xi.size());
    #pragma omp parallel for
    for (int i=0;i<xi.size();i++){
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
