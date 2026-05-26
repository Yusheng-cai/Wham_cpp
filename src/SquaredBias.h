#pragma once
#include "Bias.h"

#include <cmath>

// A simple biasing potential that acts on a single order paramter
//
// NOTES:
//   U_bias(x) is the following term:
//      U_bias(x) = phi * x^2
//
class SquaredBias: public Bias
{
    public:
        SquaredBias(const ParameterPack& pack);

        virtual Real calculate(const std::vector<Real>& x) override;
        virtual std::vector<Real> calculateForce(const std::vector<Real>& x) override;

    private:
        std::vector<Real> phi_;
        int dimension_;
};
