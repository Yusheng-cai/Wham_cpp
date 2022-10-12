#include "SimpleBias.h"

namespace BiasRegistry
{
    registry<SimpleBias> registerSimpleBias("simplebias");
} 


SimpleBias::SimpleBias(const ParameterPack& pack)
:Bias(pack)
{
    // read in dimension
    pack.ReadNumber("dimension", ParameterPack::KeyType::Required, dimension_);
    kappa_.resize(dimension_, 0.0);
    xstar_.resize(dimension_,0.0);
    phi_.resize(dimension_,0.0);

    // Read in kappa and phi
    bool kapparead = pack.ReadVectorNumber("kappa", ParameterPack::KeyType::Optional, kappa_);
    bool phiread   = pack.ReadVectorNumber("phi", ParameterPack::KeyType::Optional, phi_);

    hkappa_.resize(dimension_);
    ASSERT((kappa_.size() == dimension_), "The size of the input kappa does not match dimension");
    ASSERT((phi_.size() == dimension_), "The size of the phi input does not match dimension.");

    for (int i=0;i<dimension_;i++){
        hkappa_[i] = 0.5*kappa_[i];
    }

    if ( kapparead){
        pack.ReadVectorNumber("xstar", ParameterPack::KeyType::Required, xstar_);
        ASSERT((xstar_.size() == dimension_), "The dimension of xstar must match that of the kappa.");
    }
}

SimpleBias::Real SimpleBias::calculate(const std::vector<Real>& x)
{
    int size = x.size();
    ASSERT((size >= dimension_), "The dimension of the bias=" << xstar_.size() << " must be smaller than \
    the input data size = " << x.size());

    Real energy_ = 0.0;
    for (int i=0;i<dimension_;i++){
        energy_ += hkappa_[i]*std::pow(x[i] - xstar_[i],2.0);
        energy_ += phi_[i] * x[i];
    }

    return energy_;
}

std::vector<SimpleBias::Real> SimpleBias::calculateForce(const std::vector<Real>& x)
{
    int size = x.size();
    ASSERT((size >= dimension_), "The dimension of the bias=" << xstar_.size() << " must be smaller not \
    the input data size = " << x.size());

    std::vector<Real> force;
    force.resize(x.size(),0.0);

    for (int i=0;i<dimension_;i++)
    {
        force[i] = -2.0 * hkappa_[i] * (x[i] - xstar_[i]) - phi_[i]; 
    }

    return force;
}