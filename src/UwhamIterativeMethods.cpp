#include "UwhamIterativeMethods.h"

namespace UwhamCalculationStrategyRegistry
{
    registry<UwhamIterativeMethods> registerIterative("iterative");
};

void UwhamIterativeMethods::calculate()
{
    int Nsim = BUki_.getNR();
    int Ndata= BUki_.getNC();

    fnr_.resize(Nsim);
    fsc_.resize(Nsim);
    std::fill(fnr_.begin(), fnr_.end(), 0.0);
    std::fill(fsc_.begin(), fsc_.end(), 0.0);

    Eigen::VectorXd gradVec;
    std::vector<Real> grad = Gradient(BUki_, fk_, N_);
    gradVec = Eigen::Map<Eigen::VectorXd>(grad.data(), Nsim);

    bool converged = false;

    while (! converged)
    {
        // Find the hessian 
        Matrix<Real> hess = Hessian(BUki_, fk_, N_);

        Eigen::MatrixXd hessMat = Eigen::Map<Eigen::MatrixXd>(hess.data(), Nsim, Nsim);
        Eigen::VectorXd Hinvg = hessMat.bdcSvd(ComputeThinU | ComputeThinV).solve(gradVec);

        for (int i=0;i<Nsim;i++)
        {
            fnr_[i] = fnr_[i] - Hinvg[i];
        }

        // Normalized by the last one
        for (int i=0;i<Nsim;i++)
        {
            fnr_[i] = fnr_[i] - fnr_[Nsim-1];
        }

        // Find the gradient for nr
        auto NRgradient = Gradient(fnr_, BUki_, N_);


    }
    

}

std::vector<UwhamIterativeMethods::Real> UwhamIterativeMethods::Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
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

Matrix<UwhamIterativeMethods::Real> UwhamIterativeMethods::Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }
    int Nsim = BUki_.getNR();
    int Ndata = BUki_.getNC();

    Matrix<Real> Hessian(Nsim, Nsim);

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