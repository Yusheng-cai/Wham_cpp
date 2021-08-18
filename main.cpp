#include "src/TimeSeries.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "src/Wham.h"

#include <vector>
#include <iostream>
#include <string>
#include <memory>

int main(int argc, char** argv)
{
    using Real = CommonTypes::Real;
    using Whamptr = std::unique_ptr<Wham>;

    std::string fname(argv[1]);

    InputParser ip;
    ParameterPack pack;
    ip.ParseFile(fname, pack);

    auto WhamPack = pack.findParamPack("wham", ParameterPack::KeyType::Required);

    std::string WhamType;
    WhamPack->ReadString("type", ParameterPack::KeyType::Required, WhamType);
    Whamptr w = Whamptr(WhamRegistry::Factory::instance().create(WhamType, pack));
 
    return 0;
}