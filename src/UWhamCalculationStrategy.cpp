#include "UwhamCalculationStrategy.h"

UWhamCalculationStrategy::UWhamCalculationStrategy(UwhamStrategyInput& input)
:BUki_(input.BUki_), N_(input.N)
{
    input.pack.ReadNumber("printevery", ParameterPack::KeyType::Optional,print_every_);

    lnwji_.resize(BUki_.getNC());

    input.pack.ReadString("name", ParameterPack::KeyType::Required, name_);
};