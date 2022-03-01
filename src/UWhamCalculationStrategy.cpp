#include "UwhamCalculationStrategy.h"

UWhamCalculationStrategy::UWhamCalculationStrategy(UwhamStrategyInput& input)
:BUki_(input.BUki_), N_(input.N), fk_(input.fk)
{
    input.pack.ReadNumber("printevery", ParameterPack::KeyType::Optional,print_every_);

    lnwji_.resize(BUki_.getNC());
}