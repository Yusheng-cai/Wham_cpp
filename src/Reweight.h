#pragma once 

#include "tools/Assert.h"
#include "Wham.h"
#include "tools/OutputFunction.h"

#include <vector>
#include <array>
#include <string>
#include <memory>


struct ReweightInput
{
    Wham* wham_;
};

class Reweight
{
    public:
        using outputptr = std::unique_ptr<Output>;
        Reweight();

        virtual void calculate() = 0;
    private:
        outputptr output_;
};