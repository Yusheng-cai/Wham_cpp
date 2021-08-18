#include "Bias.h"

Bias::Bias(const ParameterPack& pack)
{
    pack.ReadNumber("temperature", ParameterPack::KeyType::Optional, temperature_);
    beta_ = 1.0/(Constants::R*temperature_);
}