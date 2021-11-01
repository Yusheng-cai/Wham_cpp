#pragma once
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/Constants.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <array>


// base class for all the biases
class Bias
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        Bias(const ParameterPack& pack);
        virtual ~Bias(){};

        virtual Real calculate(const std::vector<Real>& x) = 0;
        virtual std::vector<Real> calculateForce(const std::vector<Real>& x)=0;
        Real getBeta() const {return beta_;}

        const std::vector<Real>& getXstar() const {return xstar_;}

    protected:    
        Real temperature_ = 298.15;
        Real beta_;
        std::vector<Real> xstar_;
};

namespace BiasRegistry
{
    using Key = std::string;
    using Base= Bias;

    using Factory = GenericFactory<Base, Key, const ParameterPack&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const ParameterPack&>;
};