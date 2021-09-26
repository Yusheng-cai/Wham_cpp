#pragma once
#include "tools/GenericFactory.h"
#include "tools/CommonTypes.h"
#include "Array.h"
#include "Wham.h"
#include "Eigen/Dense"
#include "Eigen/Core"

#include <vector>
#include <string>
#include <array>

struct BwhamStrategyInput
{
    using Real = CommonTypes::Real;

    Matrix<Real>& BWil_;
    std::vector<Real>& N;

    // how much data are in each bin
    std::vector<Real>& Ml;
    ParameterPack& pack;
};

class BWhamCalculationStrategy
{
    public:
        using Real = CommonTypes::Real;

        BWhamCalculationStrategy(BwhamStrategyInput& input);
        virtual ~BWhamCalculationStrategy(){};

        virtual void calculate() = 0;
        const std::vector<Real>& getfk() const {return fk_;}
        const std::vector<Real>& getlnpl() const {return lnpl_;}

    protected:
        Matrix<Real>& BWil_;
        std::vector<Real>& N_;
        std::vector<Real>& Ml_;
        std::vector<Real> fk_;
        std::vector<Real> lnpl_;

        int Nsim_;
        int Nbins_;
};

namespace BwhamCalculationStrategyRegistry
{
    using Key = std::string;
    using Base= BWhamCalculationStrategy; 

    using Factory = GenericFactory<Base, Key, BwhamStrategyInput&>;

    template <typename D>
    using registry = RegisterInFactory<Base, D, Key, BwhamStrategyInput&>;
};