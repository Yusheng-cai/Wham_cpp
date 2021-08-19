#include "UwhamLBFGS.h"

namespace UwhamCalculationStrategyRegistry
{
    registry<UwhamLBFGS> registerLBFGS("LBFGS");
};

UwhamLBFGS::UwhamLBFGS(UwhamStrategyInput& input)
:UWhamCalculationStrategy(input)
{
    UwhamNLLInput in = {BUki_, N_};
    NLLeq_ = NLLptr(new UwhamNLL(in));
}

void UwhamLBFGS::calculate()
{
    LBFGSpp::LBFGSParam<Real> Param;
    Param.epsilon = 1e-6;
    Param.max_iterations = 100;

    LBFGSpp::LBFGSSolver<Real> solver(Param);

    Eigen::VectorXd fk = Eigen::VectorXd::Zero(BUki_.getNR());
    
    Real fx;

    solver.minimize(*NLLeq_,fk, fx);

    for (int i=0;i<fk.size();i++)
    {
        fk_[i] = fk[i] - fk[BUki_.getNR()-1];
    }

    lnwji_ = NLLeq_ -> calculatelnWi(BUki_, fk_, N_);
}

UwhamNLL::UwhamNLL(UwhamNLLInput& input)
:N_(input.N_), BUki_(input.BUki)
{
    fk_.resize(BUki_.getNR());
    N_fraction_.resize(BUki_.getNR());
    Ntot_ = 0;

    for (int i=0;i<BUki_.getNR();i++)
    {
        Ntot_ += N_[i]; 
    }

    for (int i=0;i<BUki_.getNR();i++)
    {
        N_fraction_[i] = N_[i]/Ntot_;
    }
}


UwhamNLL::Real UwhamNLL::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
    int Nsim = BUki_.getNR();
    int Ndata= BUki_.getNC();

    ASSERT((x.size() == Nsim), "The dimension of fk does not match that of the number of simulation.");

    for (int i=0;i<x.size();i++)
    {
        fk_[i] = x[i];
    }

    Real firstPart = 0.0;

    for (int i=0;i<Nsim;i++)
    {
        firstPart += N_[i] * x[i];
    }

    firstPart /= Ntot_;

    Real secondPart = 0.0;

    #pragma omp parallel for reduction(+:secondPart)
    for (int i=0;i<Ndata;i++)
    {
        std::vector<Real> temp(Nsim);
        for (int j=0;j<Nsim;j++)
        {
            temp[j] = x[j] - BUki_(j,i);
        }

        secondPart += WhamTools::LogSumExp(temp,N_fraction_);
    }

    secondPart /= Ntot_;

    Gradient(BUki_, fk_, N_, grad);

    return -firstPart + secondPart;
}

std::vector<UwhamNLL::Real> UwhamNLL::calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
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

void UwhamNLL::Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N, Eigen::VectorXd& grad)
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

    grad = Eigen::Map<Eigen::VectorXd>(gradient.data(), Nsim);
}

