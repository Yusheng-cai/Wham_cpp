#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "UwhamCalculationStrategy.h"

#include <vector>
#include <array>
#include <memory>
#include <iostream>

class Uwham:public Wham
{
    public:
        using Biasptr = std::unique_ptr<Bias>;
        using stratptr= std::unique_ptr<UWhamCalculationStrategy>;

        Uwham(const ParameterPack& pack);
        virtual ~Uwham(){};
        void initializeStrat(const ParameterPack* pack);
        void initializeBins(const std::vector<const ParameterPack*>& BinPacks);

        virtual void calculate();
        virtual void printOutput() override;

    private:
        Matrix<Real> BUki_;
        std::vector<std::vector<Real>> xi_;

        std::vector<Biasptr> Biases_;

        int Ntot_;

        stratptr strat_;

        std::string NormalizationFileOutput_;
        std::ofstream NormalizationFileofs_;

        int dimension_;

        std::vector<Bin> Bins_;
};