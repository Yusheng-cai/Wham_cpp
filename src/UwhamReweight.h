#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "Bias.h"
#include "Reweight.h"

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <map>
#include <functional>

class Uwham;

struct UwhamReweightInputPack
{
    Uwham& uwham;
    ParameterPack& pack;
};

class UwhamReweight
{
    public:
        using Real = CommonTypes::Real;
        using Biasptr = std::unique_ptr<Bias>;
        using outputfunc = std::function<void(std::string)>;

        // inputted pack in the whamPack
        UwhamReweight(UwhamReweightInputPack& pack);

        // inputs are the lnwji weights and the xi points
        void calculate();

        void printAverages(std::string name);
        void printFE(std::string name);

        void printOutput();

        void registerOutputFunc(std::string name, outputfunc func);
        outputfunc& getOutputByName(std::string name);

    private:
        std::vector<Biasptr> Vectorbias_;
        int numBias_;
        int numBins_;

        // input parameters
        Uwham& wham_;
        ParameterPack& pack_;

        // The ones vector used for calculating each of the normalizing factors 
        std::vector<Real> ones_;

        // Free energy 
        std::vector<std::vector<Real>> FE_;

        // The output files
        std::map<std::string, outputfunc> MapNameToOutput_;
        std::vector<std::string> VectorOutputs_;
        std::vector<std::string> VectorOutputFiles_;

        // The averages of each set of data under new potential
        std::vector<std::vector<Real>> averages_;

        // dimension of data 
        int dimension_;
};