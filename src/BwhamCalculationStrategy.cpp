#include "BwhamCalculationStrategy.h"

BWhamCalculationStrategy::BWhamCalculationStrategy(BwhamStrategyInput& input)
:BWil_(input.BWil_), N_(input.N), Ml_(input.Ml)
{
    Nsim_ = BWil_.getNR();
    Nbins_ = BWil_.getNC();

    fk_.resize(N_.size(),0.0);
    lnpl_.resize(Nbins_, 0.0);
}
    