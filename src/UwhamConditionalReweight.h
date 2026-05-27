#pragma once

#include "Reweight.h"
#include "Uwham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <array>
#include <string>
#include <vector>

class UwhamConditionalReweight : public Reweight
{
public:
    using Real = CommonTypes::Real;
    UwhamConditionalReweight(const ReweightInput& input);

    virtual void calculate() override;

    void printConditionalAverage(std::string name);

private:
    Uwham* uwhamptr_ = nullptr;

    // the axis at which we want to perform average conditioned on
    int axis_;

    // dimension of the data
    int dimension_;

    std::vector<std::vector<Real>> ConditionalIndex_;
    std::vector<std::vector<std::vector<Real>>> ConditionalAverage_;
};