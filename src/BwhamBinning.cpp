#include "BwhamBinning.h"

#include "tools/Assert.h"

namespace
{
using BinSelection = std::vector<int>;

int totalBinsFor(const std::vector<const Bin*>& bins)
{
    int totalBins = 1;
    for (const Bin* bin : bins)
    {
        totalBins *= bin->getNumbins();
    }
    return totalBins;
}

void appendGridPoint(
    const std::vector<const Bin*>& bins,
    int dataDimension,
    const BinSelection& selection,
    BwhamBinning::BinGrid& grid)
{
    std::vector<BwhamBinning::Real> center(dataDimension, 0.0);
    std::vector<int> binIndex(dataDimension, 0);

    for (std::size_t i=0; i<bins.size(); i++)
    {
        int dimension = bins[i]->getDimension() - 1;
        ASSERT((dimension >= 0 && dimension < dataDimension), "BWHAM bin dimension is outside the data dimension.");

        center[dimension] = bins[i]->getLocationOfBin(selection[i]);
        binIndex[dimension] = selection[i];
    }

    auto inserted = grid.indexToFlat.insert(std::make_pair(binIndex, static_cast<int>(grid.centers.size())));
    ASSERT((inserted.second), "There is duplicate BWHAM bin index.");
    grid.centers.push_back(center);
}

void buildGridRecursive(
    const std::vector<const Bin*>& bins,
    int dataDimension,
    std::size_t binNumber,
    BinSelection& selection,
    BwhamBinning::BinGrid& grid)
{
    if (binNumber == bins.size())
    {
        appendGridPoint(bins, dataDimension, selection, grid);
        return;
    }

    for (int i=0; i<bins[binNumber]->getNumbins(); i++)
    {
        selection[binNumber] = i;
        buildGridRecursive(bins, dataDimension, binNumber + 1, selection, grid);
    }
}
}

BwhamBinning::BinGrid BwhamBinning::buildBinGrid(const std::vector<const Bin*>& bins, int dataDimension)
{
    ASSERT((!bins.empty()), "BWHAM requires at least one bin definition.");
    ASSERT((static_cast<int>(bins.size()) == dataDimension), "The number of BWHAM bin definitions must match the data dimension.");
    ASSERT((dataDimension <= 2), "BWHAM currently supports at most two dimensions.");

    BinGrid grid;
    grid.totalBins = totalBinsFor(bins);
    grid.centers.reserve(grid.totalBins);

    BinSelection selection(bins.size(), 0);
    buildGridRecursive(bins, dataDimension, 0, selection, grid);

    return grid;
}

bool BwhamBinning::findBinIndexForSample(
    const std::vector<const Bin*>& bins,
    const std::vector<Real>& sample,
    std::vector<int>& binIndex)
{
    ASSERT((sample.size() == bins.size()), "BWHAM sample dimension must match the number of bin definitions.");

    binIndex.assign(bins.size(), 0);

    for (const Bin* bin : bins)
    {
        int dimension = bin->getDimension() - 1;
        ASSERT((dimension >= 0 && dimension < static_cast<int>(sample.size())), "BWHAM bin dimension is outside the sample dimension.");

        Real value = sample[dimension];
        if (!bin->isInRange(value))
        {
            return false;
        }

        binIndex[dimension] = bin->findBin(value);
    }

    return true;
}
