#include "UwhamCalculationStrategy.h"
#include "parallel/OpenMP_buffer.h"

#include <vector>
#include <array>

class UwhamAdaptiveMethods: public UWhamCalculationStrategy
{
    public:
        UwhamAdaptiveMethods(UwhamStrategyInput& input);
        virtual ~UwhamAdaptiveMethods(){};

        virtual void calculate();

        Matrix<Real> Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
        std::vector<Real> Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
        std::vector<Real> calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);
        Real calculateError(const std::vector<Real>& fi, const std::vector<Real>& fi_prev);

    private:
        std::vector<Real> fnr_;
        std::vector<Real> fsc_;

        // tolerance is default to 1e-7
        Real tolerance_ = 1e-7;
};