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
        const std::vector<Real>& getNorm() const {return norms_;}

        virtual void calculate(std::vector<Real>& fk) = 0;
        std::string getName() {return name_;}

    protected:
        Matrix<Real>& BUki_;
        std::vector<Real>& N_;
        std::vector<Real> fk_;
        std::vector<Real> lnwji_;
        std::vector<Real> norms_;

        // print frequency
        int print_every_=-1;
        std::string name_;
};

namespace UwhamCalculationStrategyRegistry
{
    using Key = std::string;
    using Base= UWhamCalculationStrategy; 

    using Factory = GenericFactory<Base, Key, UwhamStrategyInput&>;

    template <typename D>
    using registry = RegisterInFactory<Base, D, Key, UwhamStrategyInput&>;
};