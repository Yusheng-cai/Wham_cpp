#pragma once
#include "BWham.h"
#include "BwhamCalculationStrategy.h"
#include "tools/CommonTypes.h"
#include "LBFGS/LBFGS.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

struct BwhamNLLInput
{
    using Real = CommonTypes::Real;
    Matrix<Real>& BWil;
    std::vector<Real>& N_; 
    std::vector<Real>& Ml_;
};

class BwhamNLL
{
    public:
        using Real = CommonTypes::Real;

        BwhamNLL(BwhamNLLInput& input);

        Real operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad);

    private:
        Matrix<Real>& BWil_;
        std::vector<Real>& N_;
        std::vector<Real> N_fraction_;
        std::vector<Real>& Ml_;
        Real Ntot_ = 0.0;
        std::vector<Real> fk_;
};

class BwhamLBFGS : public BWhamCalculationStrategy
{
    public:
        using NLLptr = std::unique_ptr<BwhamNLL>;

        BwhamLBFGS(BwhamStrategyInput& input);

        virtual void calculate();
    
    private:
        NLLptr NLLeq_;

        Real epsilon_ = 1e-6;

        int max_iterations_ = 1e5;
};