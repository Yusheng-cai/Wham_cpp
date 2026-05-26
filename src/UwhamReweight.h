#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "Bias.h"
#include "Uwham.h"
#include "Reweight.h"

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <map>
#include <functional>

class UwhamReweight : public Reweight
{
    public:
        using Real = CommonTypes::Real;
        using Biasptr = std::unique_ptr<Bias>;
        using outputfunc = std::function<void(std::string)>;

        // inputted pack in the whamPack
        UwhamReweight(const ReweightInput& input);

        // inputs are the lnwji weights and the xi points
        void calculate();

        void printReweightAverages(std::string name);

        void printFreeEnergys(std::string name);

    private:
        // input parameters
        Uwham* Uwham_;

        // The averages of each set of data under new potential
        std::vector<std::vector<Real>> averages_;

        // dimension of data 
        int dimension_;

        std::vector<std::map<std::vector<int>, Real>> FreeEnergys_;
};