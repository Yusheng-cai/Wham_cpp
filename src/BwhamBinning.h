#pragma once

#include "Bin.h"
#include "tools/CommonTypes.h"

#include <map>
#include <vector>

namespace BwhamBinning
{
    using Real = CommonTypes::Real;

    struct BinGrid
    {
        int totalBins = 0;
        std::vector<std::vector<Real>> centers;
        std::map<std::vector<int>, int> indexToFlat;
    };

    BinGrid buildBinGrid(const std::vector<const Bin*>& bins, int dataDimension);
    bool findBinIndexForSample(
        const std::vector<const Bin*>& bins,
        const std::vector<Real>& sample,
        std::vector<int>& binIndex);
}
