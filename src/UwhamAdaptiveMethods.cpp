#include "UwhamAdaptiveMethods.h"

namespace UwhamCalculationStrategyRegistry
{
    registry<UwhamAdaptiveMethods> registerAdaptive("adaptive");
};

UwhamAdaptiveMethods::UwhamAdaptiveMethods(UwhamStrategyInput& input)
:UWhamCalculationStrategy(input)
{
    input.pack.ReadNumber("tolerance",ParameterPack::KeyType::Optional, tolerance_);
}

void UwhamAdaptiveMethods::calculate()
{
    int Nsim = BUki_.getNR();
    int Ndata= BUki_.getNC();

    fnr_.resize(Nsim);
    fsc_.resize(Nsim);
    std::vector<Real> ones(Ndata);
    std::fill(ones.begin(), ones.end(), 1.0);
    std::fill(fnr_.begin(), fnr_.end(), 0.0);
    std::fill(fsc_.begin(), fsc_.end(), 0.0);

    Eigen::VectorXd gradVec;


    bool converged = false;

    Real err = 0.0;

    while ( ! converged)
    {
        std::vector<Real> grad = Gradient(BUki_, fk_, N_);
        gradVec = Eigen::Map<Eigen::VectorXd>(grad.data(), Nsim);

        // Find the hessian 
        Matrix<Real> hess = Hessian(BUki_, fk_, N_);

        Eigen::MatrixXd hessMat = Eigen::Map<Eigen::MatrixXd>(hess.data(), Nsim, Nsim);
        Eigen::VectorXd Hinvg = hessMat.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(gradVec);

        for (int i=0;i<Nsim;i++)
        {
            fnr_[i] = fk_[i] - Hinvg[i];
        }

        // Normalized by the last one
        for (int i=0;i<Nsim;i++)
        {
            fnr_[i] = fnr_[i] - fnr_[Nsim-1];
        }

        // Find the gradient for nr
        auto NRgradient = Gradient(BUki_,fnr_, N_);
        Real normNR = WhamTools::NormVector(NRgradient);

        std::vector<Real> lnwji = calculatelnWi(BUki_, fk_, N_);
        for (int i=0;i<Nsim;i++)
        {
            std::vector<Real> column;
            column.resize(Ndata);
            #pragma omp parallel for
            for (int j=0;j<Ndata;j++)
            {
                column[j] = lnwji[j] - BUki_(i,j);
            }
            fsc_[i] = -1.0*WhamTools::LogSumExp(column, ones);
        }

        for (int i=0;i<Nsim;i++)
        {
            fsc_[i] = fsc_[i] - fsc_[Nsim-1];
        }

        auto SCgradient = Gradient(BUki_, fsc_, N_);  
        Real normSC = WhamTools::NormVector(SCgradient);

        if (normSC < normNR)
        {
            err = calculateError(fsc_, fk_);
            fk_.assign(fsc_.begin(), fsc_.end());
        }
        else
        {
            err = calculateError(fnr_, fk_);
            fk_.assign(fnr_.begin(), fnr_.end());
        }

        if (err < tolerance_)
        {
            converged = true;
        }
    }

    lnwji_ = calculatelnWi(BUki_, fk_, N_);
}

UwhamAdaptiveMethods::Real UwhamAdaptiveMethods::calculateError(const std::vector<Real>& fi, const std::vector<Real>& fi_prev)
{
    ASSERT((fi.size() == fi_prev.size()), "The sizes don't match in calculating error.");
    std::vector<Real> absdiff_;
    std::vector<Real> absprev_;
    absdiff_.resize(fi.size());
    absprev_.resize(fi.size());

    for (int i=0;i<absdiff_.size();i++)
    {
        absdiff_[i] = std::abs(fi[i] - fi_prev[i]);
    }

    for (int i=0;i<absprev_.size();i++)
    {
        absprev_[i] = std::abs(fi_prev[i]);
    }

    auto it = std::max_element(absdiff_.begin(), absdiff_.end()-1);
    auto itprev = std::max_element(absprev_.begin(), absprev_.end()-1);

    return (*it)/(*itprev);
}

std::vector<UwhamAdaptiveMethods::Real> UwhamAdaptiveMethods::calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{ 
    int Nsim = BUki_.getNR();
    int Ndata = BUki_.getNC();

    std::vector<Real> lnwji;
    lnwji.resize(Ndata);

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<Ndata;i++)
        {
            std::vector<Real> column;
            column.resize(Nsim);
            for (int j=0;j<Nsim;j++)
            {
                column[j] = fk[j]-1.0*BUki_(j,i);
            }

            Real val = WhamTools::LogSumExp(column, N);
            lnwji[i] = -1.0*val;
        }
    }

    return lnwji;
}

std::vector<UwhamAdaptiveMethods::Real> UwhamAdaptiveMethods::Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }

    int Nsim = BUki_.getNR();
    int Ndata= BUki_.getNC();

    std::vector<Real> ones_(Ndata);
    std::fill(ones_.begin(), ones_.end(), 1.0);

    std::vector<Real> gradient(Nsim);
    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N); 

    std::vector<std::vector<Real>> lnpki(Nsim, std::vector<Real>(Ndata));

    #pragma omp parallel for collapse(2)
    for (int k=0;k<Nsim;k++)
    {
        for (int j=0;j<Ndata;j++)
        {
            lnpki[k][j] = fk[k] - BUki_(k,j) + lnwji[j];
        }
    }

    std::vector<Real> lnpk(Nsim);

    for (int k=0;k<Nsim;k++)
    {
        Real val = WhamTools::LogSumExp(lnpki[k], ones_);
        lnpk[k] = val;
    }

    for (int k=0;k<Nsim;k++)
    {
        gradient[k] = -1.0/Ntot*(N[k] - N[k] * std::exp(lnpk[k]));
    }

    return gradient;
}


Matrix<UwhamAdaptiveMethods::Real> UwhamAdaptiveMethods::Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }
    int Nsim = BUki_.getNR();
    int Ndata = BUki_.getNC();

    Matrix<Real> Hessian(Nsim, Nsim);

    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N); 
    

    std::vector<std::vector<Real>> pki(Nsim, std::vector<Real>(Ndata));

    #pragma omp parallel for collapse(2)
    for (int k=0;k<Nsim;k++)
    {
        for (int j=0;j<Ndata;j++)
        {
            Real lnpki;
            lnpki = fk[k] - BUki_(k,j) + lnwji[j];
            pki[k][j] = std::exp(lnpki);
        }
    }

    for (int i=0;i<Nsim;i++)
    {
        for (int j=0;j<Nsim;j++)
        {
            if (i == j)
            {
                Real sum = 0.0;
                Real sum_sq = 0.0;
                #pragma omp parallel for reduction(+:sum,sum_sq) 
                for (int k=0;k<Ndata;k++)
                {
                    sum += pki[i][k];
                    sum_sq += pki[i][k] * pki[i][k];
                }

                Hessian(i,j) = -1.0/Ntot*(-N[i]*sum + N[i]*N[i]*sum_sq);
            }
            else
            {
                Real sum = 0.0;
                #pragma omp parallel for reduction(+:sum)
                for (int k=0;k<Ndata;k++)
                {
                    sum += pki[i][k] * pki[j][k];
                }

                Hessian(i,j) = -1.0/Ntot*(sum*N[i]*N[j]);
            }
        }
    }

    return Hessian;
}