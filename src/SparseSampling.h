#pragma once 

#include "Wham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>

class SparseSampling : public Wham
{
    public:
        SparseSampling(WhamInput& input);

        // make the assumption that input bias and timeseries has the same order
        virtual void calculate() override;

    private:
        std::vector<std::vector<Real>> means_;
        std::vector<std::vector<Real>> std_;
        std::vector<Real> energy_;
        std::vector<Real> force_;
};