#pragma once 

#include "tools/OutputFunction.h"
#include "TimeSeries.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

struct TSInput
{
    using tsptr = std::shared_ptr<TimeSeries>;
    ParameterPack& pack_;
    std::vector<tsptr>& vectorTS_;
};

class TSoperation
{
    public:
        using tsptr = std::shared_ptr<TimeSeries>;
        using OutputFuncPtr = std::unique_ptr<Output>;

        TSoperation(const TSInput& input);

        virtual void calculate() = 0;
        virtual void print();
    
    protected:
        std::vector<tsptr>& VectorTS_;

        OutputFuncPtr outputs_;

        ParameterPack& pack_;

        std::vector<std::string> VectorOutputNames_;
        std::vector<std::string> VectorOutputFileNames_;
};

namespace timeseriesOP
{
    using Base = TSoperation;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,const TSInput&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const TSInput&>;
}