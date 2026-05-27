#pragma once

#include "Bin.h"
#include "TSoperation.h"
#include "parallel/OpenMP.h"
#include "parallel/OpenMP_buffer.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

class HistogramOP : public TSoperation
{
public:
    using Binptr = std::unique_ptr<Bin>;
    using Real = CommonTypes::Real;
    HistogramOP(const TSInput& input);

    virtual void calculate() override;
    void printHistogram(std::string name);
    void printTotalHistogram(std::string name);

private:
    std::vector<Binptr> Bins_;
    int binDimension_;
    std::map<std::vector<int>, int> histogramMap_;

    std::vector<std::vector<std::vector<Real>>> histogram_;
    std::map<std::vector<int>, int> MapIndexToHistogram_;
};