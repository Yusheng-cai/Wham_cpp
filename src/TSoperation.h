#pragma once 
#include "tools/OutputFunction.h"
#include "TimeSeries.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

struct TSInput
{
    std::vector<TimeSeries>& vectorTS_;
    ParameterPack& pack_;
};

class TSoperation
{
    public:
        using OutputFuncPtr = std::unique_ptr<Output>;

        TSoperation(const TSInput& input);

        virtual void calculate();
        virtual void print();
    
    private:
        std::vector<TimeSeries>& VectorTS_;

        OutputFuncPtr outputs_;

        ParameterPack& pack_;
};