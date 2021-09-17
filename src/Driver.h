#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"

#include <memory>
#include <vector>

class Driver
{
    public:
        using Whamptr = std::unique_ptr<Wham>;
        using tsptr   = std::shared_ptr<TimeSeries>;

        Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

        void InitializeWham(const ParameterPack& pack);
        void calculate();
        void printOutput();

    private:
        std::vector<Whamptr> VectorWhamCalc_;
        std::vector<tsptr> VectorTimeSeries_;
};