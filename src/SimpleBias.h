#pragma once
#include "Bias.h"

#include <cmath>

// A simple biasing potential that acts on a single order paramter
// 
// NOTES: 
//   U_bias(x) is the sum of the following terms:
//
//     U_harmonic(x)   = 0.5*kappa*(x - x_star)^2 
//     U_linear        = phi*x + constant
class SimpleBias: public Bias
{
    public:
        SimpleBias(const ParameterPack& pack);

        virtual Real calculate(const std::vector<Real>& x) override;

    private:
        std::vector<Real> kappa_;
        std::vector<Real> hkappa_;
        std::vector<Real> phi_;
        std::vector<Real> xstar_;
        int dimension_;
};