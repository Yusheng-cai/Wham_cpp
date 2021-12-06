#pragma once

#include "Reweight.h"
#include "tools/Assert.h"
#include "Uwham.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>

class UwhamConditionReweight : public Reweight
{
    public:
        using Real = CommonTypes::Real;
        UwhamConditionReweight(const ReweightInput& input);

        virtual void calculate() override;

        void printAverage(std::string name);

    private:
        Uwham* uwhamptr_=nullptr;

        int axis_;
        int axisavg_;

        std::vector<std::vector<Real>> axisdata_;
        std::vector<std::vector<Real>> peraxisaverage_;
};