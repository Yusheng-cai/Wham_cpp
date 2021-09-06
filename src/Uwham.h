#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "UwhamCalculationStrategy.h"

#include <vector>
#include <iomanip>
#include <map>
#include <array>
#include <memory>
#include <iostream>
#include <chrono>

class Uwham:public Wham
{
    public:
        using Biasptr = std::unique_ptr<Bias>;
        using stratptr= std::unique_ptr<UWhamCalculationStrategy>;

        Uwham(const WhamInput& pack);
        virtual ~Uwham(){};
        void initializeStrat(const ParameterPack* pack);
        void initializeBins(const std::vector<const ParameterPack*>& BinPacks);
        void OpenFile(std::ofstream& ofs, std::string& name);

        virtual void calculate() override;
        virtual void printOutput() override;

    private:
        Matrix<Real> BUki_;
        std::vector<std::vector<Real>> xi_;

        std::vector<Biasptr> Biases_;

        int Ntot_;

        stratptr strat_;

        std::string NormalizationFileOutput_;
        std::ofstream NormalizationFileofs_;

        std::string pjiFileOutput_;
        std::ofstream pjiFileofs_;

        std::string lnwjiOutput_;
        std::ofstream lnwjiFileofs_;

        int dimension_;

        std::vector<Bin> Bins_;

        // The lnwji that falls within each of the bins
        std::map<std::vector<int>, std::vector<Real>> MapBinIndexToVectorlnwji_;
        std::map<std::vector<int>, Real> MapBinIndexToWji_;

        // precision of the output
        int precision_=3;
};