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

        // Calculate the error between 2 iterations
        Real calculateError(const std::vector<Real>& fi, const std::vector<Real>& fi_prev);

        // Calculate the Newton Raphson step
        void NewtonRaphsonStep(std::vector<Real>& fnr, std::vector<Real>& gradientNR);

        // calculate the self-consistent step
        void SelfConsistentStep(std::vector<Real>& fsc, std::vector<Real>& gradientSC);


    private:
        std::vector<Real> fnr_;
        std::vector<Real> fsc_;

        // tolerance is default to 1e-7
        Real tolerance_ = 1e-7;

        // print frequency
        int print_every_=-1;
};