#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"
#include "TSoperation.h"
#include "Reweight.h"

#include <memory>
#include <map>
#include <vector>

class Driver
{
    public:
        using Whamptr = std::unique_ptr<Wham>;
        using tsptr   = std::shared_ptr<TimeSeries>;
        using TSopptr = std::unique_ptr<TSoperation>;
        using reweightptr = std::unique_ptr<Reweight>;

        Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

        void InitializeWham();
        void InitializeTSoperation();
        void InitializeReweight();
        void calculate();
        void finishCalculate();
        void printOutput();

    private:
        std::vector<Whamptr> VectorWhamCalc_;
        std::vector<tsptr> VectorTimeSeries_;
        std::vector<TSopptr> VectorTimeSeriesOP_;
        std::vector<reweightptr> ReweightPtr_;
        ParameterPack& pack_;

        // map the name of the wham calculation to the location of it in the vector
        std::map<std::string, int> MapNameOfWhamToLoc_;
};