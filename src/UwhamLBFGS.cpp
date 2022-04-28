#include "UwhamLBFGS.h"

namespace UwhamCalculationStrategyRegistry
{
    registry<UwhamLBFGS> registerLBFGS("LBFGS");
};

UwhamLBFGS::UwhamLBFGS(UwhamStrategyInput& input)
:UWhamCalculationStrategy(input)
{
    input.pack.ReadNumber("epsilon", ParameterPack::KeyType::Optional, epsilon_);
    input.pack.ReadNumber("max_iterations", ParameterPack::KeyType::Optional, max_iterations_);
    input.pack.ReadNumber("epsilon_rel", ParameterPack::KeyType::Optional, epsilon_rel_);

    UwhamNLLInput in = {BUki_, N_};
    NLLeq_ = NLLptr(new UwhamNLL(in));
}

void UwhamLBFGS::calculate(std::vector<Real>& fk)
{
    LBFGSpp::LBFGSParam<Real> Param;
    Param.epsilon = epsilon_;
    Param.max_iterations = max_iterations_;
    Param.epsilon_rel = epsilon_rel_;

    LBFGSpp::LBFGSSolver<Real> solver(Param);

    // copy the data -> fk
    fk_ = fk;
    Eigen::VectorXd fk_E = Eigen::Map<Eigen::VectorXd>(fk.data(), BUki_.getNR());
    Real fx;

    int numiterations = solver.minimize(*NLLeq_,fk_E, fx, print_every_);
    norms_ = NLLeq_->getNorms();

    for (int i=0;i<fk_E.size();i++)
    {
        fk_[i] = fk_E[i];
    }

    // subtract the minimum fk 
    Real normalize = fk_[0];
    for (int i=0;i<fk_.size();i++)
    {
        fk_[i] = fk_[i] - normalize;
    }

    lnwji_ = WhamTools::calculatelnWi(BUki_, fk_, N_);

    // need to reweight lnwji
    std::vector<Real> ones(lnwji_.size(),1);
    Real f = -1.0*WhamTools::LogSumExpOMP(lnwji_, ones);

    #pragma omp parallel for 
    for (int i=0;i<lnwji_.size();i++)
    {
        lnwji_[i] = f + lnwji_[i];
    }

    for (int i=0;i<fk_.size();i++)
    {
        fk_[i] = fk_[i] - f;
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

    std::vector<Real> fk(x.size(),0.0);

    for (int i=0;i<x.size();i++)
    {
        fk[i] = x[i] - x[0];
    }

    Real value = WhamTools::Uwham_NLL_equation(fk, BUki_, N_);
    auto gradient = WhamTools::Gradient(BUki_, fk, N_);
    grad = Eigen::Map<Eigen::VectorXd>(gradient.data(), Nsim); 

    derives_.push_back(grad);
    norms_.push_back(grad.norm());

    return value;
}
