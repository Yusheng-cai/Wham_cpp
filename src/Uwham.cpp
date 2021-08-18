#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const ParameterPack& input)
:Wham(input)
{
    N_.resize(VectorTimeSeries_.size());
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xji_.insert(xji_.end(),VectorTimeSeries_[i].begin(), VectorTimeSeries_[i].end());
        N_[i] = VectorTimeSeries_[i].getSize(); 

        Ntot_ += N_[i];
    }

    auto biases = input.findParamPacks("bias", ParameterPack::KeyType::Required);
    ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

    for (int i=0;i<biases.size();i++)
    {
        std::string biastype = "simplebias";
        biases[i] -> ReadString("type", ParameterPack::KeyType::Optional, biastype);
        Biasptr b = Biasptr(BiasRegistry::Factory::instance().create(biastype, *biases[i]));

        Biases_.push_back(std::move(b));
    }

    BUji_.resize(biases.size(), xji_.size());
    for (int i=0;i<biases.size();i++)
    {
        for (int j=0;j<xji_.size();j++)
        {
            Real val = Biases_[i]->calculate(xji_[j]); 
            BUji_(i,j) = val;
        }
    }
}

void Uwham::calculate()
{

}

Matrix<Uwham::Real> Uwham::Hessian(const Matrix<Real>& BUji, const std::vector<Real>& fi, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }

    Matrix<Real> Hessian(BUji_.getNR(), BUji_.getNR());

    std::vector<Real> lnwji;
    lnwji.resize(BUji_.getNC());

    for (int i=0;i<lnwji.size();i++)
    {
        std::vector<Real> column;
        column.resize(BUji_.getNR());
        for (int j=0;j<BUji_.getNR();j++)
        {
            column[j] = fi[j]-1.0*BUji_(j,i);
        }

        Real val = WhamTools::LogSumExp(column, N);
        lnwji[i] = val;
    }
}