#include "UwhamCalculationStrategy.h"

UWhamCalculationStrategy::UWhamCalculationStrategy(UwhamStrategyInput& input)
:BUki_(input.BUki_), N_(input.N), fk_(input.fk)
{
    lnwji_.resize(BUki_.getNC());
}