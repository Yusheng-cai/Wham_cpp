#pragma once 

#include "tools/Assert.h"
#include "tools/OutputFunction.h"
#include "tools/InputParser.h"
#include "Bias.h"
#include "tools/GenericFactory.h"
#include "Wham.h"

#include <vector>
#include <array>
#include <string>
#include <memory>


struct ReweightInput
{
    Wham* wham_;
    ParameterPack& pack_;
};

class Reweight
{
    public:
        using outputfunc = std::function<void(std::string)>;
        using outputptr = std::unique_ptr<Output>;
        using Biasptr    = std::unique_ptr<Bias>;
        Reweight(const ReweightInput& input);

        void printOutput();
        void checkOutputValidity();

        virtual void calculate() = 0;
    protected:
        outputptr output_;
        ParameterPack& pack_;

        // keep a record of the WHAM object 
        Wham* wham_;

        std::vector<std::string> outputNames_;
        std::vector<std::string> outputfileNames_;

        int numBias_;
        std::vector<Biasptr> Vectorbias_;
};

namespace ReweightRegistry
{
    using Base = Reweight;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,const ReweightInput&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const ReweightInput&>;
};