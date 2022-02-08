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

void UwhamLBFGS::calculate()
{
    LBFGSpp::LBFGSParam<Real> Param;
    Param.epsilon = epsilon_;
    Param.max_iterations = max_iterations_;
    Param.epsilon_rel = epsilon_rel_;

    LBFGSpp::LBFGSSolver<Real> solver(Param);

    // copy the data -> fk
    Eigen::VectorXd fk = Eigen::VectorXd::Zero(BUki_.getNR());
    for (int i=0;i<BUki_.getNR();i++)
    {
        fk[i] = fk_[i];
    }
    
    Real fx;

    int numiterations = solver.minimize(*NLLeq_,fk, fx);
    std::cout << "It took " << numiterations << " iterations to converge." << "\n";

    for (int i=0;i<fk.size();i++)
    {
        //fk_[i] = fk[i] - fk[BUki_.getNR()-1];
        fk_[i] = fk[i];
    }
    Real min = *(std::min_element(fk_.begin(), fk_.end()));
    for (int i=0;i<fk_.size();i++)
    {
        fk_[i] = fk_[i] - min;
    }

    std::cout << "Function value = " << fx << "\n";

    lnwji_ = WhamTools::calculatelnWi(BUki_, fk_, N_);

    // need to reweigth lnwji
    std::vector<Real> ones(lnwji_.size(),1);
    Real f = -1.0*WhamTools::LogSumExpOMP(lnwji_, ones);

    for (int i=0;i<fk.size();i++)
    {
        //fk_[i] = fk[i] - fk[BUki_.getNR()-1];
        fk_[i] = fk_[i] - f;
    }

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
        fk_[i] = x[i] - x[0];
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
