#include "UwhamLBFGS.h"

UwhamLBFGS::UwhamLBFGS(UwhamStrategyInput& input)
:UWhamCalculationStrategy(input)
{
    input.pack.ReadNumber("epsilon", ParameterPack::KeyType::Optional, epsilon_);
    input.pack.ReadNumber("max_iterations", ParameterPack::KeyType::Optional, max_iterations_);
    input.pack.ReadNumber("epsilon_rel", ParameterPack::KeyType::Optional, epsilon_rel_);

    UwhamNLLInput in = {BUki_, N_};
    NLLeq_ = NLLptr(new UwhamNLL(in));
}

UwhamStrategyResult UwhamLBFGS::calculate(const std::vector<Real>& fk)
{
    LBFGSpp::LBFGSParam<Real> Param;
    Param.epsilon = epsilon_;
    Param.max_iterations = max_iterations_;
    Param.epsilon_rel = epsilon_rel_;

    LBFGSpp::LBFGSSolver<Real> solver(Param);

    // copy the data -> fk
    std::vector<Real> fk_result = fk;
    Eigen::VectorXd fk_E = Eigen::Map<const Eigen::VectorXd>(fk.data(), BUki_.getNR());
    Real fx;

    // do the solver minimization
    solver.minimize(*NLLeq_,fk_E, fx, print_every_);
    std::vector<Real> norms = NLLeq_->getNorms();

    for (int i=0;i<fk_E.size();i++){
        fk_result[i] = fk_E[i];
    }

    // subtract the minimum fk 
    Real normalize = fk_result[0];
    for (int i=0;i<fk_result.size();i++){
        fk_result[i] = fk_result[i] - normalize;
    }

    std::vector<Real> lnwji = WhamTools::calculatelnWi(BUki_, fk_result, N_);

    // need to reweight lnwji
    std::vector<Real> ones(lnwji.size(),1);
    Real f = -1.0*WhamTools::LogSumExpOMP(lnwji, ones);

    #pragma omp parallel for 
    for (int i=0;i<lnwji.size();i++){
        lnwji[i] = f + lnwji[i];
    }

    // normalize fk --> -log(Qi/Q0)
    fk_result = fk_result - f;

    return {fk_result, lnwji, norms};
}

UwhamNLL::UwhamNLL(UwhamNLLInput& input)
:BUki_(input.BUki), N_(input.N_)
{
    fk_.resize(BUki_.getNR());

    N_fraction_.resize(BUki_.getNR());
    Ntot_ = 0;

    for (int i=0;i<BUki_.getNR();i++){
        Ntot_ += N_[i]; 
    }

    // Ni / Ntot
    N_fraction_ = N_ / Ntot_;
}


UwhamNLL::Real UwhamNLL::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad){
    int Nsim = BUki_.getNR();

    ASSERT((x.size() == Nsim), "The dimension of fk does not match that of the number of simulation.");

    std::vector<Real> fk(x.size(),0.0);

    for (int i=0;i<x.size();i++){fk[i] = x[i] - x[0];}

    Real value = WhamTools::Uwham_NLL_equation(fk, BUki_, N_);
    auto gradient = WhamTools::Gradient(BUki_, fk, N_);
    grad = Eigen::Map<Eigen::VectorXd>(gradient.data(), Nsim); 

    derives_.push_back(grad);
    norms_.push_back(grad.norm());

    return value;
}
