#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const ParameterPack& input)
:Wham(input)
{
    N_.resize(VectorTimeSeries_.size());
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTimeSeries_[i].begin(), VectorTimeSeries_[i].end());
        N_[i] = VectorTimeSeries_[i].getSize(); 

        Ntot_ += N_[i];
    }

    auto biases = input.findParamPacks("bias", ParameterPack::KeyType::Required);
    ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

    for (int i=0;i<biases.size();i++)
    {
        std::string biastype = "simplebias";
        biases[i] -> ReadString("type", ParameterPack::KeyType::Optional, biastype);
        Biasptr b = Biasptr(BiasRegistry::Factory::instance().create(biastype, *biases[i]));

        Biases_.push_back(std::move(b));
    }

    BUki_.resize(biases.size(), xi_.size());
    for (int i=0;i<biases.size();i++)
    {
        for (int j=0;j<xi_.size();j++)
        {
            Real val = Biases_[i]->calculate(xi_[j]); 
            BUki_(i,j) = Biases_[i]->getBeta()*val;
        }
    }

    // Now read in what type of calculation are you letting Uwham do
    auto whamPack = input.findParamPack("wham", ParameterPack::KeyType::Required);
    initializeStrat(whamPack);
}

void Uwham::initializeStrat(const ParameterPack* whampack)
{
    std::string strat;
    whampack->ReadString("strategy", ParameterPack::KeyType::Required, strat);

    UwhamStrategyInput input = {BUki_, N_};
    strat_ = stratptr(UwhamCalculationStrategyRegistry::Factory::instance().create(strat, input));
}

void Uwham::calculate()
{
    strat_ -> calculate();
}
