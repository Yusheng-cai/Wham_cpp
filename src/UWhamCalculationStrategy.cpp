#include "UwhamCalculationStrategy.h"

UWhamCalculationStrategy::UWhamCalculationStrategy(UwhamStrategyInput& input)
:BUki_(input.BUki_), N_(input.N), fk_(input.fk)
{
    input.pack.ReadNumber("printevery", ParameterPack::KeyType::Optional,print_every_);

    lnwji_.resize(BUki_.getNC());

    input.pack.ReadString("name", ParameterPack::KeyType::Required, name_);
};

void UWhamCalculationStrategy::setFk(std::vector<Real>& fk)
{
    ASSERT((fk.size() == BUki_.getNR()), "The set fk does not match the size");
    fk_.clear();
    fk_.insert(fk_.end(), fk.begin(), fk.end());
}