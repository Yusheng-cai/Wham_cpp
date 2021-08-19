#pragma once
#include "Eigen/Dense"
#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"
#include "tools/Constants.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <array>
#include <string>
#include <chrono>

struct WhamInput
{
    ParameterPack& pack_;
    std::vector<TimeSeries>& VectorTimeSeries_;
};

class Wham
{
    public:
        using Real = CommonTypes::Real;

        Wham(const WhamInput& input);
        virtual ~Wham(){};

        virtual void calculate() = 0;
        virtual void printOutput() {};
    
    protected:
        std::vector<TimeSeries>& VectorTimeSeries_;

        std::vector<Real> N_;
};

namespace WhamRegistry
{
    using Base = Wham;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,const WhamInput&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const WhamInput&>;
};

namespace WhamTools
{
    using Real = CommonTypes::Real;

    // Pass in vector is calculated as Log(\sum(N*exp(vector)))
    Real LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N);

    // find the norm of a vector
    Real NormVector(const std::vector<Real>& vector);

    // find the hessian matrix of the UWham NLL equation
    Matrix<Real> Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    // find the gradient vector of the UWham NLL equation
    std::vector<Real> Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    // find the lnWi in UWham
    std::vector<Real> calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
};