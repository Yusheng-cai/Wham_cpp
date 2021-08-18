#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"

#include <memory>
#include <vector>

class Driver
{
    public:
        using Whamptr = std::unique_ptr<Wham>;

        Driver(const ParameterPack& pack);

        void InitializeWham(const ParameterPack& pack);

    private:
        Whamptr whamCalc_;
};