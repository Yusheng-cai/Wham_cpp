#include "Wham.h"

Wham::Wham(const ParameterPack& pack)
{
    // find all the instances of the timeseries block
    auto TsPacks = pack.findParamPacks("Timeseries", ParameterPack::KeyType::Required);

    for (int i= 0 ;i < TsPacks.size();i++)
    {
        VectorTimeSeries_.push_back(TimeSeries(*TsPacks[i]));
    }

    // currently let's just assume that we are dealing with one wham pack at a time
    auto whamPack = pack.findParamPack("wham", ParameterPack::KeyType::Required);

    // what information should be contained in the wham pack?
    
    // temperature is assumed to be 298.15 if not specified
    whamPack -> ReadNumber("temperature", ParameterPack::KeyType::Optional, temperature_);
    beta_ = 1.0/(Constants::R*temperature_);
}