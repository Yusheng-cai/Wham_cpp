#include "UwhamCalculationStrategy.h"
#include "parallel/OpenMP_buffer.h"

#include <vector>
#include <array>

class UwhamIterativeMethods: public UWhamCalculationStrategy
{
    public:
        UwhamIterativeMethods(UwhamStrategyInput& input):UWhamCalculationStrategy(input){};
        virtual ~UwhamIterativeMethods(){};

        virtual void calculate();

        Matrix<Real> Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
        std::vector<Real> Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    private:
        std::vector<Real> fnr_;
        std::vector<Real> fsc_;
};