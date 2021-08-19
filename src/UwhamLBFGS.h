#include "Uwham.h"
#include "tools/CommonTypes.h"
#include "LBFGS/LBFGS.h"
#include "Eigen/Dense"
#include "Eigen/Core"

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

        void Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N, Eigen::VectorXd& grad);
        std::vector<Real> calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
    
    private:
        Matrix<Real>& BUki_;
        std::vector<Real>& N_;
        std::vector<Real> N_fraction_;
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
};