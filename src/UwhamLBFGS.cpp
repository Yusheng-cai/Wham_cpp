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

    lnwji_ = WhamTools::calculatelnWi(BUki_, fk_, N_);

    // need to reweigth lnwji
    std::vector<Real> ones(lnwji_.size(),1);
    Real f = -1.0*WhamTools::LogSumExp(lnwji_, ones);

    #pragma omp parallel for 
    for (int i=0;i<lnwji_.size();i++)
    {
        lnwji_[i] = f + lnwji_[i];
    }
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

    auto gradient = WhamTools::Gradient(BUki_, fk_, N_);

    grad = Eigen::Map<Eigen::VectorXd>(gradient.data(), Nsim); 

    return -firstPart + secondPart;
}
