#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"
#include "TSoperation.h"

#include <memory>
#include <vector>

class Driver
{
    public:
        using Whamptr = std::unique_ptr<Wham>;
        using tsptr   = std::shared_ptr<TimeSeries>;
        using TSopptr = std::unique_ptr<TSoperation>;

        Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

        void InitializeWham(const ParameterPack& pack);
        void InitializeTSoperation(const ParameterPack& pack);
        void calculate();
        void finishCalculate();
        void printOutput();

    private:
        std::vector<Whamptr> VectorWhamCalc_;
        std::vector<tsptr> VectorTimeSeries_;
        std::vector<TSopptr> VectorTimeSeriesOP_;
};