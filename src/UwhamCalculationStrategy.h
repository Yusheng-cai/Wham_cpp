#pragma once

#include "Array.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Wham.h"
#include "tools/CommonTypes.h"

#include <array>
#include <string>
#include <vector>

struct UwhamStrategyInput
{
    using Real = CommonTypes::Real;

    Matrix<Real>& BUki_;
    std::vector<Real>& N;
    ParameterPack& pack;
};

struct UwhamStrategyResult
{
    using Real = CommonTypes::Real;

    std::vector<Real> fk;
    std::vector<Real> lnwji;
    std::vector<Real> norms;
};

class UWhamCalculationStrategy
{
public:
    using Real = CommonTypes::Real;

    UWhamCalculationStrategy(UwhamStrategyInput& input);
    virtual ~UWhamCalculationStrategy() {};

    virtual UwhamStrategyResult calculate(const std::vector<Real>& fk) = 0;
    std::string getName() { return name_; }

protected:
    Matrix<Real>& BUki_;
    std::vector<Real>& N_;

    // print frequency
    int print_every_ = -1;
    std::string name_;
};
