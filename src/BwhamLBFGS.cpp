#include "BwhamLBFGS.h"

namespace BwhamCalculationStrategyRegistry
{
    registry<BwhamLBFGS> registerLBFGS("LBFGS");
};

BwhamLBFGS::BwhamLBFGS(BwhamStrategyInput& input)
:BWhamCalculationStrategy(input)
{
    input.pack.ReadNumber("epsilon", ParameterPack::KeyType::Optional, epsilon_);
    input.pack.ReadNumber("max_iterations", ParameterPack::KeyType::Optional, max_iterations_);

    BwhamNLLInput in = {BWil_, N_, Ml_};
    NLLeq_ = NLLptr(new BwhamNLL(in));
}

void BwhamLBFGS::calculate()
{
    LBFGSpp::LBFGSParam<Real> Param;
    Param.epsilon = epsilon_;
    Param.max_iterations = max_iterations_;

    LBFGSpp::LBFGSSolver<Real> solver(Param);

    Eigen::VectorXd fk = Eigen::VectorXd::Zero(BWil_.getNR());
    
    Real fx;

    solver.minimize(*NLLeq_,fk, fx);

    for (int i=0;i<fk.size();i++)
    {
        fk_[i] = fk[i] - fk[BWil_.getNR()-1];
        std::cout << "fk " << i << " = " << fk_[i] << std::endl;
    }

    // obtain the lnpl
    lnpl_ = WhamTools::calculatelnpl(BWil_, Ml_, N_, fk_);

    // lnwji_ = WhamTools::calculatelnWi(BWil_, fk_, N_);

    // need to reweigth lnwji
    //std::vector<Real> ones(lnwji_.size(),1);
    //Real f = -1.0*WhamTools::LogSumExpOMP(lnwji_, ones);

    // #pragma omp parallel for 
    // for (int i=0;i<lnwji_.size();i++)
    // {
    //     lnwji_[i] = f + lnwji_[i];
    // }
}

BwhamNLL::BwhamNLL(BwhamNLLInput& input)
:N_(input.N_), BWil_(input.BWil) , Ml_(input.Ml_)
{
    fk_.resize(BWil_.getNR(),0.0);
    N_fraction_.resize(BWil_.getNR());
    Ntot_ = 0;

    for (int i=0;i<BWil_.getNR();i++)
    {
        Ntot_ += N_[i]; 
    }

    for (int i=0;i<BWil_.getNR();i++)
    {
        N_fraction_[i] = N_[i]/Ntot_;
    }
}


BwhamNLL::Real BwhamNLL::operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
    int Nsim = BWil_.getNR();
    int Ndata= BWil_.getNC();

    ASSERT((x.size() == Nsim), "The dimension of fk does not match that of the number of simulation.");

    for (int i=0;i<x.size();i++)
    {
        fk_[i] = x[i];
        fk_[i]  -= x[Nsim-1];
    }

    Real firstPart = 0.0;

    for (int i=0;i<Nsim;i++)
    {
        firstPart -= N_[i] * fk_[i];
    }

    Real secondPart = 0.0;

    for (int i=0;i<Ndata;i++)
    {
        std::vector<Real> temp(Nsim);
        for (int j=0;j<Nsim;j++)
        {
            temp[j] = fk_[j] - BWil_(j,i);
        }

        secondPart += -Ml_[i] * std::log(Ml_[i]) + Ml_[i] * WhamTools::LogSumExp(temp,N_);
    }

    std::cout << "Total = " << firstPart + secondPart << std::endl;

    auto gradient = WhamTools::BGradient(BWil_, Ml_, N_, fk_);

    grad = Eigen::Map<Eigen::VectorXd>(gradient.data(), Nsim); 

    // This step is especially important when the bin does not cover all the data points
    Real sum_grad_=0.0;
    for (int i=0;i<grad.size();i++)
    {
        sum_grad_ -= grad[i];
    }

    grad[Nsim-1] +=  sum_grad_;

    return firstPart + secondPart;
}