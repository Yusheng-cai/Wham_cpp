#pragma once

#include "TSoperation.h"

#include <vector>
#include <string>
#include <array>

class CombineData : public TSoperation
{
    public:
        CombineData(const TSInput& input);

        virtual void calculate() override {};
    private:
};