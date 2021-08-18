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
};

class UWhamCalculationStrategy
{
    public:
        using Real = CommonTypes::Real;

        UWhamCalculationStrategy(UwhamStrategyInput& input);
        virtual ~UWhamCalculationStrategy(){};

        virtual void calculate() = 0;

    protected:
        Matrix<Real>& BUki_;
        std::vector<Real>& N_;
        std::vector<Real> fk_;
};

namespace UwhamCalculationStrategyRegistry
{
    using Key = std::string;
    using Base= UWhamCalculationStrategy; 

    using Factory = GenericFactory<Base, Key, UwhamStrategyInput&>;

    template <typename D>
    using registry = RegisterInFactory<Base, D, Key, UwhamStrategyInput&>;
};