#include "tools/GenericFactory.h"
#include "tools/CommonTypes.h"
#include "Array.h"
#include "Wham.h"
#include "Eigen/Dense"
#include "Eigen/Core"

#include <vector>
#include <string>
#include <array>

struct UwhamStrategyInput
{
    using Real = CommonTypes::Real;

    Matrix<Real>& BUki_;
    std::vector<Real>& N;
    ParameterPack& pack;
    std::vector<Real> fk;
};

class UWhamCalculationStrategy
{
    public:
        using Real = CommonTypes::Real;

        UWhamCalculationStrategy(UwhamStrategyInput& input);
        virtual ~UWhamCalculationStrategy(){};

        const std::vector<Real>& getFk_() const {return fk_;}
        const std::vector<Real>& getlnwji_() const {return lnwji_;}

        virtual void calculate() = 0;

    protected:
        Matrix<Real>& BUki_;
        std::vector<Real>& N_;
        std::vector<Real> fk_;
        std::vector<Real> lnwji_;
};

namespace UwhamCalculationStrategyRegistry
{
    using Key = std::string;
    using Base= UWhamCalculationStrategy; 

    using Factory = GenericFactory<Base, Key, UwhamStrategyInput&>;

    template <typename D>
    using registry = RegisterInFactory<Base, D, Key, UwhamStrategyInput&>;
};