#include "SimpleBias.h"

namespace BiasRegistry
{
    registry<SimpleBias> registerSimpleBias("simplebias");
} 


SimpleBias::SimpleBias(const ParameterPack& pack)
:Bias(pack)
{
    bool kapparead = pack.ReadVectorNumber("kappa", ParameterPack::KeyType::Optional, kappa_);
    dimension_ = kappa_.size();
    hkappa_.resize(dimension_);

    for (int i=0;i<dimension_;i++)
    {
        hkappa_[i] = 0.5*kappa_[i];
    }

    if ( kapparead)
    {
        pack.ReadVectorNumber("xstar", ParameterPack::KeyType::Required, xstar_);
        ASSERT((xstar_.size() == dimension_), "The dimension of xstar must match that of the kappa.");
    }
}

SimpleBias::Real SimpleBias::calculate(const std::vector<Real>& x)
{
    int size = x.size();
    ASSERT((x.size() == dimension_), "The dimension of the bias=" << xstar_.size() << " and does not \
    the input data size = " << x.size());

    Real energy_ = 0.0;
    for (int i=0;i<size;i++)
    {
        energy_ += hkappa_[i]*std::pow(x[i] - xstar_[i],2.0);
    }

    return energy_;
}