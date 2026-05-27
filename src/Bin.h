#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <array>
#include <cmath>
#include <vector>

class Bin
{
public:
    using Real = CommonTypes::Real;
    using Range = CommonTypes::Real2;

    Bin(const ParameterPack& pack);

    int findBin(Real data) const;
    bool isInRange(Real data) const;

    Range getRange() const { return range_; }
    int getNumbins() const { return numbins_; }
    Real getStep() const { return step_; }
    int getDimension() const { return dimension_; }
    Real getLocationOfBin(int binNum) const;

private:
    Range range_;
    int numbins_;
    int dimension_ = 1;
    Real step_;
};
