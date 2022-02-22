#include "UwhamAdaptiveMethods.h"

namespace UwhamCalculationStrategyRegistry
{
    registry<UwhamAdaptiveMethods> registerAdaptive("adaptive");
};

UwhamAdaptiveMethods::UwhamAdaptiveMethods(UwhamStrategyInput& input)
:UWhamCalculationStrategy(input)
{
    input.pack.ReadNumber("tolerance",ParameterPack::KeyType::Optional, tolerance_);
    input.pack.ReadNumber("printevery", ParameterPack::KeyType::Optional,print_every_);
}

void UwhamAdaptiveMethods::NewtonRaphsonStep(std::vector<Real>& fnr, std::vector<Real>& gradientNR)
{
    int Nsim = BUki_.getNR();

    std::vector<Real> grad = WhamTools::Gradient(BUki_, fk_, N_);
    Eigen::VectorXd gradVec = Eigen::Map<Eigen::VectorXd>(grad.data(), Nsim);

    fnr.clear();
    fnr.resize(Nsim,0.0);

    // Find the hessian 
    Matrix<Real> hess = WhamTools::Hessian(BUki_, fk_, N_);
    Eigen::MatrixXd hessMat = Eigen::Map<Eigen::MatrixXd>(hess.data(), Nsim, Nsim);
    Eigen::VectorXd Hinvg = (hessMat.transpose() * hessMat).ldlt().solve(hessMat.transpose() * gradVec);

    for (int i=0;i<Nsim;i++)
    {
        fnr[i] = fk_[i] - Hinvg[i];
    }

    // normalize by the first value
    Real normalizedNR = fnr[0];
    for (int i=0;i<Nsim;i++)
    {
        fnr[i] = fnr[i] - normalizedNR;
    }

    gradientNR = WhamTools::Gradient(BUki_,fnr, N_);
}

void UwhamAdaptiveMethods::SelfConsistentStep(std::vector<Real>& fsc, std::vector<Real>& gradientSC)
{
    int Nsim  = BUki_.getNR();
    int Ndata = BUki_.getNC();
    std::vector<Real> ones(Ndata , 1.0);

    fsc.clear();
    fsc.resize(Nsim, 0.0);

    std::vector<Real> lnwji = WhamTools::calculatelnWi(BUki_, fk_, N_);
    for (int i=0;i<Nsim;i++)
    {
        std::vector<Real> column;
        column.resize(Ndata);
        #pragma omp parallel for
        for (int j=0;j<Ndata;j++)
        {
            column[j] = lnwji[j] - BUki_(i,j);
        }
        fsc[i] = -1.0*WhamTools::LogSumExpOMP(column, ones);
    }

    Real SC_normalized = fsc[0];
    for (int i=0;i<Nsim;i++)
    {
        fsc[i] = fsc[i] - SC_normalized;
    }

    gradientSC = WhamTools::Gradient(BUki_, fsc, N_);  
}

void UwhamAdaptiveMethods::calculate()
{
    int Nsim = BUki_.getNR();
    int Ndata= BUki_.getNC();

    std::vector<Real> NRgradient(Nsim, 0.0);
    std::vector<Real> SCgradient(Nsim, 0.0);

    fnr_.resize(Nsim,0.0);
    fsc_.resize(Nsim,0.0);
    std::vector<Real> ones(Ndata, 1.0);
 
    Eigen::VectorXd gradVec;

    bool converged = false;

    Real err = 0.0;
    int step = 1;

    while ( ! converged)
    {
        // update Newton Raphson using current guess of fk
        NewtonRaphsonStep(fnr_, NRgradient);
        Real normNR = WhamTools::NormVector(NRgradient);

        // update self consistent using current guess of fk 
        SelfConsistentStep(fsc_, SCgradient); 
        Real normSC = WhamTools::NormVector(SCgradient);

        if (normSC < normNR)
        {
            err = calculateError(fsc_, fk_);
            fk_.assign(fsc_.begin(), fsc_.end());
            norms_.push_back(std::sqrt(normSC));
        }
        else
        {
            err = calculateError(fnr_, fk_);
            fk_.assign(fnr_.begin(), fnr_.end());
            norms_.push_back(std::sqrt(normNR));
        }

        if ((print_every_ != -1) && (step % print_every_==0)) 
        {
            std::cout << "self consistent norm is " << normSC << " NR norm is " << normNR << "\n";

            if (normSC < normNR)
            {
                std::cout << "Self consistent iteration is chosen because of lower norm." << "\n";
            }
            else
            {
                std::cout << "Newton Raphson is chosen because of lower norm." << "\n";
            }
            std::cout << "The error at step " << step << " is " << err << std::endl;
        }

        if (err < tolerance_)
        {
            converged = true;
        }

        step++;
    }

    lnwji_ = WhamTools::calculatelnWi(BUki_, fk_, N_);

    // need to reweight lnwji
    Real f = -1.0*WhamTools::LogSumExp(lnwji_, ones);

    #pragma omp parallel for 
    for (int i=0;i<lnwji_.size();i++)
    {
        lnwji_[i] = f + lnwji_[i];
    }

    // calculates the NLL equation
    Real NLL_val = WhamTools::Uwham_NLL_equation(fk_, BUki_, N_);

    std::cout << "NLL value = " << NLL_val << "\n";
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