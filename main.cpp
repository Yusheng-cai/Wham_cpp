#include "src/TimeSeries.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <vector>
#include <iostream>

int main(int argc, char** argv)
{
    using Real = CommonTypes::Real;

    std::string fname(argv[1]);

    InputParser ip;
    ParameterPack pack;
    ip.ParseFile(fname, pack);

    auto tspack = pack.findParamPack("timeseries", ParameterPack::KeyType::Required);

    TimeSeries ts(*tspack);

    return 0;
}