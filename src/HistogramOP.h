#pragma once 

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "TSoperation.h"
#include "Bin.h"

#include <vector>
#include <array>
#include <memory>
#include <string>
#include <map>

class HistogramOP : public TSoperation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Real   = CommonTypes::Real;
        HistogramOP(const TSInput& input);

        virtual void calculate() override;
        void printHistogram(std::string name);

    private:
        std::vector<Binptr> Bins_;
        int binDimension_;
        std::vector<std::vector<Real>> xi_;
        std::map<std::vector<int>, int> histogramMap_;
};