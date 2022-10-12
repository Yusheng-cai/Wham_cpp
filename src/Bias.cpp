#include "Bias.h"

Bias::Bias(const ParameterPack& pack){
    pack.ReadNumber("temperature", ParameterPack::KeyType::Optional, temperature_);
    // in kJ/mol
    beta_ = 1000.0/(Constants::R*temperature_);
}