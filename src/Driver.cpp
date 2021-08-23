#include "Driver.h"

Driver::Driver(const ParameterPack& pack, const CommandLineArguments& cmd)
{
    // read in the absolute path of the timeseries input, if not read then timeseries will assume it is the current dir
    std::string abspath;
    bool read = cmd.readString("abspath", CommandLineArguments::Keys::Optional,abspath);

    // find all the instances of the timeseries block
    auto TsPacks = pack.findParamPacks("timeseries", ParameterPack::KeyType::Required);

    for (int i= 0 ;i < TsPacks.size();i++)
    {
        TimeSeriesInputPack input = { const_cast<ParameterPack&>(*TsPacks[i]), abspath};
        VectorTimeSeries_.push_back(TimeSeries(input));
    }
    InitializeWham(pack);
}

void Driver::InitializeWham(const ParameterPack& pack)
{
    auto whamPack = pack.findParamPack("wham", ParameterPack::KeyType::Required);

    std::string whamType;
    whamPack->ReadString("type", ParameterPack::KeyType::Required,whamType);

    WhamInput input = {const_cast<ParameterPack&>(pack), VectorTimeSeries_};
    whamCalc_ = Whamptr(WhamRegistry::Factory::instance().create(whamType, input));
}

void Driver::calculate()
{
    whamCalc_ -> calculate();
}

void Driver::printOutput()
{
    whamCalc_ -> printOutput();
}