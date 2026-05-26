#include "SquaredBias.h"

namespace BiasRegistry
{
    registry<SquaredBias> registerSquaredBias("squaredbias");
}

SquaredBias::SquaredBias(const ParameterPack& pack)
:Bias(pack)
{
    // read in dimension
    pack.ReadNumber("dimension", ParameterPack::KeyType::Required, dimension_);
    xstar_.resize(dimension_,0.0);
    phi_.resize(dimension_,0.0);

    // Read in kappa and phi
    pack.ReadVectorNumber("phi", ParameterPack::KeyType::Optional, phi_);
    ASSERT((static_cast<int>(phi_.size()) == dimension_), "The size of the phi input does not match dimension.");
}

SquaredBias::Real SquaredBias::calculate(const std::vector<Real>& x)
{
    int size = x.size();
    ASSERT((size >= dimension_), "The dimension of the bias=" << xstar_.size() << " must be smaller than \
    the input data size = " << x.size());

    Real energy_ = 0.0;
    for (int i=0;i<dimension_;i++){
        energy_ += phi_[i] * x[i] * x[i];
    }

    return energy_;
}

std::vector<SquaredBias::Real> SquaredBias::calculateForce(const std::vector<Real>& x)
{
    int size = x.size();
    ASSERT((size >= dimension_), "The dimension of the bias=" << xstar_.size() << " must be smaller not \
    the input data size = " << x.size());

    std::vector<Real> force;
    force.resize(x.size(),0.0);

    for (int i=0;i<dimension_;i++){
        force[i] = -2.0 * phi_[i] * x[i];
    }

    return force;
}
