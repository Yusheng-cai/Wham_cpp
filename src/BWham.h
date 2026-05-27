#pragma once
#include "Array.h"
#include "Bin.h"
#include "BwhamCalculationStrategy.h"
#include "Wham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <iomanip>
#include <map>
#include <memory>
#include <vector>

class Bwham : public Wham
{
public:
    using BinPtr = std::unique_ptr<Bin>;
    using StrategyPtr = std::unique_ptr<BWhamCalculationStrategy>;

    Bwham(const WhamInput& input);

    virtual void calculate() override;
    virtual std::string type() override { return "Bwham"; }

    void printlnpl(std::string name);

private:
    void initializeBinnedGrid();
    void countDataPerBin();
    void initializeWil();
    void initializeStrategy();

    std::vector<BinPtr> bins_;
    std::vector<std::vector<Real>> centerBins_;

    int totalBins_ = 0;

    Matrix<Real> reducedBias_;

    std::map<std::vector<int>, int> binIndexToFlat_;

    std::vector<Real> countsPerBin_;

    StrategyPtr strategy_;

    std::vector<Real> lnpl_;
};
