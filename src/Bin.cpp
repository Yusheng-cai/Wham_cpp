#include "Bin.h"

Bin::Bin(const ParameterPack& pack)
{
    pack.ReadArrayNumber("range", ParameterPack::KeyType::Required, range_);
    pack.ReadNumber("numbins", ParameterPack::KeyType::Required,numbins_);

    step_ = (range_[1] - range_[0])/numbins_;
}

int Bin::findBin(Real x)
{
    int index = std::floor((x - range_[0])/step_);

    return index;
}