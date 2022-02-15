#pragma once
#include "Uwham.h"
#include "Wham.h"
#include "tools/CommonTypes.h"
#include "LBFGS/LBFGS.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

struct UwhamNLLInput
{
    using Real = CommonTypes::Real;
    Matrix<Real>& BUki;
    std::vector<Real>& N_; 
};

class UwhamNLL
{
    public:
        using Real = CommonTypes::Real;

        UwhamNLL(UwhamNLLInput& input);

        Real operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);

        const std::vector<Eigen::VectorXd>& getDerivatives() const { return derives_;}
        const std::vector<Real>& getNorms() const {return norms_;}

    private:
        Matrix<Real>& BUki_;
        std::vector<Real>& N_;
        std::vector<Real> N_fraction_;

        // a vector that stores the derivativses
        std::vector<Eigen::VectorXd> derives_;
        std::vector<Real> norms_;
        Real Ntot_ = 0.0;
        std::vector<Real> fk_;
};

class UwhamLBFGS : public UWhamCalculationStrategy
{
    public:
        using NLLptr = std::unique_ptr<UwhamNLL>;

        UwhamLBFGS(UwhamStrategyInput& input);

        virtual void calculate();
    
    private:
        NLLptr NLLeq_;

        Real epsilon_ = 1e-6;

        Real epsilon_rel_ = 1e-6;

        int max_iterations_ = 100;
};