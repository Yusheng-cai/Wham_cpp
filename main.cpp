#include "src/TimeSeries.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "src/Driver.h"
#include "tools/Assert.h"

#include <vector>
#include <iostream>
#include <string>
#include <memory>

int main(int argc, char** argv)
{
    using Real = CommonTypes::Real;
    using Whamptr = std::unique_ptr<Wham>;

    ASSERT((argc >= 2), "The input file must be provided.");
    std::string fname(argv[1]);

    InputParser ip;
    ParameterPack pack;
    ip.ParseFile(fname, pack);

    Driver d(pack);
    d.calculate();

 
    return 0;
}