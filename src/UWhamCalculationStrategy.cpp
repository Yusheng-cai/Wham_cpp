#include "UwhamCalculationStrategy.h"

UWhamCalculationStrategy::UWhamCalculationStrategy(UwhamStrategyInput& input)
:BUki_(input.BUki_), N_(input.N)
{
    fk_.resize(N_.size());
    std::fill(fk_.begin(), fk_.end(), 1e-8);

    lnwji_.resize(BUki_.getNC());
}