#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"

#include <memory>
#include <vector>

class Driver
{
    public:
        using Whamptr = std::unique_ptr<Wham>;

        Driver(const ParameterPack& pack);

        void InitializeWham(const ParameterPack& pack);
        void calculate();
        void printOutput();

    private:
        Whamptr whamCalc_;
};