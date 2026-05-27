#pragma once
#include "Reweight.h"
#include "TSoperation.h"
#include "TimeSeries.h"
#include "Wham.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <map>
#include <memory>
#include <vector>

class Driver
{
public:
    using WhamPtr = std::unique_ptr<Wham>;
    using TimeSeriesPtr = std::shared_ptr<TimeSeries>;
    using TimeSeriesOperationPtr = std::unique_ptr<TSoperation>;
    using ReweightPtr = std::unique_ptr<Reweight>;

    Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

    void calculate();
    void finishCalculate();
    void printOutput();

private:
    void initializeWham();
    void initializeTimeSeriesOperations();
    void initializeReweight();

    std::vector<WhamPtr> whamCalculations_;
    std::vector<TimeSeriesPtr> timeSeries_;
    std::vector<TimeSeriesOperationPtr> timeSeriesOperations_;
    std::vector<ReweightPtr> reweights_;
    ParameterPack& pack_;

    std::map<std::string, int> whamNameToIndex_;
};
