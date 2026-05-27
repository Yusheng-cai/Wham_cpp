#pragma once

#include "TSoperation.h"

#include <array>
#include <string>
#include <vector>

class CombineData : public TSoperation
{
public:
    CombineData(const TSInput& input);

    virtual void calculate() override {};

private:
};