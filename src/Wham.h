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

class Wham
{
    public:
        using Real = CommonTypes::Real;

        Wham(const ParameterPack& pack);
        virtual ~Wham(){};

        virtual void calculate() = 0;
    
    protected:
        std::vector<TimeSeries> VectorTimeSeries_;

        std::vector<int> N_;
};

namespace WhamRegistry
{
    using Base = Wham;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,const ParameterPack&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const ParameterPack&>;
};

namespace WhamTools
{
    using Real = CommonTypes::Real;

    // Pass in vector is calculated as Log(\sum(N*exp(vector)))
    Real LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N);
};