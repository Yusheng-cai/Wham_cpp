#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "UwhamCalculationStrategy.h"

#include <vector>
#include <array>
#include <memory>

class Uwham:public Wham
{
    public:
        using Biasptr = std::unique_ptr<Bias>;
        using stratptr= std::unique_ptr<UWhamCalculationStrategy>;

        Uwham(const ParameterPack& pack);
        virtual ~Uwham(){};
        void initializeStrat(const ParameterPack* pack);

        virtual void calculate();

    private:
        Matrix<Real> BUki_;
        std::vector<Real> fk_;
        std::vector<std::vector<Real>> xi_;

        std::vector<Biasptr> Biases_;

        int Ntot_;

        stratptr strat_;
};