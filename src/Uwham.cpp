#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const ParameterPack& input)
:Wham(input)
{
    N_.resize(VectorTimeSeries_.size());
    std::vector<int> dimensions_(VectorTimeSeries_.size());

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTimeSeries_[i].begin(), VectorTimeSeries_[i].end());
        N_[i] = VectorTimeSeries_[i].getSize(); 
        dimensions_[i] = VectorTimeSeries_[i].getDimension();

        Ntot_ += N_[i];
    }

    for (int i=0;i<dimensions_.size()-1;i++)
    {
        ASSERT((dimensions_[i] == dimensions_[i+1]), "The dimension in the " << i << "th timeseries does not match with the " << i+1 << "th time series");
    }

    // record the dimensions of this Wham calculation
    dimension_ = dimensions_[0];

    auto biases = input.findParamPacks("bias", ParameterPack::KeyType::Required);
    ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

    for (int i=0;i<biases.size();i++)
    {
        std::string biastype = "simplebias";
        biases[i] -> ReadString("type", ParameterPack::KeyType::Optional, biastype);
        Biasptr b = Biasptr(BiasRegistry::Factory::instance().create(biastype, *biases[i]));

        Biases_.push_back(std::move(b));
    }

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

    // Now read in what type of calculation are you letting Uwham do
    auto whamPack = input.findParamPack("wham", ParameterPack::KeyType::Required);
    initializeStrat(whamPack);

    bool normRead = whamPack->ReadString("normalizationOutput", ParameterPack::KeyType::Optional, NormalizationFileOutput_);

    if (normRead)
    {
        NormalizationFileofs_.open(NormalizationFileOutput_);
        ASSERT((NormalizationFileofs_.is_open()), "The file with name " << NormalizationFileOutput_ << " is not opened.");
    }

    auto binPacks = whamPack -> findParamPacks("bins", ParameterPack::KeyType::Required);
    initializeBins(binPacks);

    ASSERT((binPacks.size() == dimension_), "The binning dimension is " << binPacks.size() << " while the dimension of the Wham is " << dimension_);
}

void Uwham::initializeBins(const std::vector<const ParameterPack*>& BinPacks)
{
    if (BinPacks.size() != 0)
    {
        for (int i=0;i<BinPacks.size();i++)
        {
            Bins_.push_back(Bin(*BinPacks[i])); 
        }
    }
}


void Uwham::initializeStrat(const ParameterPack* whampack)
{
    std::string strat;
    whampack->ReadString("strategy", ParameterPack::KeyType::Required, strat);

    UwhamStrategyInput input = {BUki_, N_, const_cast<ParameterPack&>(*whampack)};
    strat_ = stratptr(UwhamCalculationStrategyRegistry::Factory::instance().create(strat, input));
}

void Uwham::calculate()
{
    auto start = std::chrono::high_resolution_clock::now();
    strat_ -> calculate();

    auto end = std::chrono::high_resolution_clock::now();

    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Calculation took " << diff.count() << std::endl;

    const auto& lnwji = strat_ -> getlnwji_();
}

void Uwham::printOutput()
{
    if (NormalizationFileofs_.is_open())
    {
        int Nsim = BUki_.getNR();
        NormalizationFileofs_ << "# normalization constants" << "\n";
        for (int i=0;i<Nsim;i++)
        {
            NormalizationFileofs_ << strat_ -> getFk_()[i] << "\n";
        }

        NormalizationFileofs_.close();
    }

}
