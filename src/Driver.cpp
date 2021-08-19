#include "Driver.h"

Driver::Driver(const ParameterPack& pack)
{
    InitializeWham(pack);
}

void Driver::InitializeWham(const ParameterPack& pack)
{
    auto whamPack = pack.findParamPack("wham", ParameterPack::KeyType::Required);

    std::string whamType;
    whamPack->ReadString("type", ParameterPack::KeyType::Required,whamType);

    whamCalc_ = Whamptr(WhamRegistry::Factory::instance().create(whamType, pack));
}

void Driver::calculate()
{
    whamCalc_ -> calculate();
}

void Driver::printOutput()
{
    whamCalc_ -> printOutput();
}