#include "WhamTools.h"

#include "VectorOperations.h"
#include "tools/Assert.h"

#include <algorithm>
#include <cmath>

WhamTools::Real WhamTools::LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N)
{
    ASSERT((vector.size() == N.size()), "The size of the vector is not equal to the size of N.");

    // Find the max of the vector
    auto it = std::max_element(vector.begin(), vector.end());
    Real maxVal = *it;

    Real sum = 0.0;
    for (int i = 0; i < vector.size(); i++) {
        sum += N[i] * std::exp(vector[i] - maxVal);
    }

    sum = std::log(sum) + maxVal;

    return sum;
}

WhamTools::Real WhamTools::LogSumExpOMP(const std::vector<Real>& vector, const std::vector<Real>& N)
{
    ASSERT((vector.size() == N.size()), "The size of the vector is not equal to the size of N.");

    // Find the max of the vector
    auto it = std::max_element(vector.begin(), vector.end());
    Real maxVal = *it;

    Real sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < vector.size(); i++) {
        sum += N[i] * std::exp(vector[i] - maxVal);
    }

    sum = std::log(sum) + maxVal;

    return sum;
}

WhamTools::Real WhamTools::NormVector(const std::vector<Real>& vector)
{
    Real sum_ = 0.0;

    for (int i = 0; i < vector.size(); i++) {
        sum_ += vector[i] * vector[i];
    }

    return sum_;
}

std::vector<WhamTools::Real> WhamTools::calculatelnpl(const Matrix<Real>& BWil,
                                                      const std::vector<Real>& Ml,
                                                      const std::vector<Real>& N,
                                                      const std::vector<Real>& fk)
{
    int Nbins = Ml.size();
    int Nsim = BWil.getNR();

    std::vector<Real> lnpl(Nbins, 0.0);
    std::vector<Real> temp;

    for (int i = 0; i < Nbins; i++) {
        std::vector<Real> temp(Nsim, 0.0);
        for (int j = 0; j < Nsim; j++) {
            temp[j] = fk[j] - BWil(j, i);
        }

        Real res = WhamTools::LogSumExp(temp, N);
        lnpl[i] = std::log(Ml[i]) - res;
    }

    return lnpl;
}

std::vector<WhamTools::Real> WhamTools::calculatelnWi(const Matrix<Real>& BUki,
                                                      const std::vector<Real>& fk,
                                                      const std::vector<Real>& N)
{
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    std::vector<Real> lnwji;
    lnwji.resize(Ndata);

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < Ndata; i++) {
            std::vector<Real> column;
            column.resize(Nsim);
            for (int j = 0; j < Nsim; j++) {
                column[j] = fk[j] - 1.0 * BUki(j, i);
            }

            Real val = WhamTools::LogSumExp(column, N);
            lnwji[i] = -1.0 * val;
        }
    }

    return lnwji;
}
std::vector<WhamTools::Real> WhamTools::BGradient(const Matrix<Real>& BWil,
                                                  const std::vector<Real>& Ml,
                                                  const std::vector<Real>& N,
                                                  const std::vector<Real>& fk)
{
    int Nbins = Ml.size();
    int Nsim = BWil.getNR();

    Real Ntot = 0;
    for (int i = 0; i < Nsim; i++) {
        Ntot += N[i];
    }

    // calculate lnpl
    std::vector<Real> lnpl = WhamTools::calculatelnpl(BWil, Ml, N, fk);

    std::vector<Real> derivative(Nsim, 0.0);
    for (int i = 0; i < Nsim; i++) {
        std::vector<Real> temp(Nbins, 0.0);
        std::vector<Real> ones(Nbins, 1.0);
        for (int j = 0; j < Nbins; j++) {
            temp[j] = lnpl[j] - BWil(i, j);
        }

        Real res = WhamTools::LogSumExp(temp, ones);
        derivative[i] = N[i] * (std::exp(fk[i] + res) - 1);
    }

    return derivative;
}

std::vector<WhamTools::Real> WhamTools::Gradient(const Matrix<Real>& BUki,
                                                 const std::vector<Real>& fk,
                                                 const std::vector<Real>& N)
{
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    std::vector<Real> ones_(Ndata);
    std::fill(ones_.begin(), ones_.end(), 1.0);

    std::vector<Real> gradient(Nsim);
    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N);

    std::vector<std::vector<Real>> lnpki(Nsim, std::vector<Real>(Ndata));

#pragma omp parallel for collapse(2)
    for (int k = 0; k < Nsim; k++) {
        for (int j = 0; j < Ndata; j++) {
            lnpki[k][j] = fk[k] - BUki(k, j) + lnwji[j];
        }
    }

    std::vector<Real> lnpk(Nsim);

    for (int k = 0; k < Nsim; k++) {
        Real val = WhamTools::LogSumExp(lnpki[k], ones_);
        lnpk[k] = val;
    }

    for (int k = 0; k < Nsim; k++) {
        // gradient[k] = -1.0/Ntot*(N[k] - N[k] * std::exp(lnpk[k]));
        gradient[k] = -(N[k] - N[k] * std::exp(lnpk[k]));
    }

    return gradient;
}

Matrix<WhamTools::Real> WhamTools::Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk,
                                           const std::vector<Real>& N)
{
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    Matrix<Real> Hessian(Nsim, Nsim);

    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N);

    std::vector<std::vector<Real>> pki(Nsim, std::vector<Real>(Ndata));

#pragma omp parallel for collapse(2)
    for (int k = 0; k < Nsim; k++) {
        for (int j = 0; j < Ndata; j++) {
            Real lnpki;
            lnpki = fk[k] - BUki(k, j) + lnwji[j];
            pki[k][j] = std::exp(lnpki);
        }
    }

    for (int i = 0; i < Nsim; i++) {
        for (int j = 0; j < Nsim; j++) {
            if (i == j) {
                Real sum = 0.0;
                Real sum_sq = 0.0;
#pragma omp parallel for reduction(+ : sum, sum_sq)
                for (int k = 0; k < Ndata; k++) {
                    sum += pki[i][k];
                    sum_sq += pki[i][k] * pki[i][k];
                }

                Hessian(i, j) = -1.0 * (-N[i] * sum + N[i] * N[i] * sum_sq);
            } else {
                Real sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
                for (int k = 0; k < Ndata; k++) {
                    sum += pki[i][k] * pki[j][k];
                }

                Hessian(i, j) = -1.0 * (sum * N[i] * N[j]);
            }
        }
    }

    return Hessian;
}

WhamTools::Real WhamTools::CalculateBAR(const std::vector<Real>& w_F, const std::vector<Real>& w_B,
                                        Real DeltaF)
{
    Real sizeWF = w_F.size();
    Real sizeWB = w_B.size();

    Real M = std::log(sizeWF / sizeWB);

    // log f(W) = - log [1 + exp((M + W - DeltaF))]
    //           = - log ( exp[+maxarg] [exp[-maxarg] + exp[(M + W - DeltaF) - maxarg]] )
    //           = - maxarg - log(exp[-maxarg] + exp[(M + W - DeltaF) - maxarg])
    // where maxarg = max((M + W - DeltaF), 0)
    std::vector<Real> logf_F(sizeWF, 0.0);
    std::vector<Real> onesF(sizeWF, 1.0);
#pragma omp parallel for
    for (int i = 0; i < (int)sizeWF; i++) {
        Real val = M + w_F[i] - DeltaF;
        Real maxarg = std::max(val, 0.0);

        logf_F[i] = -maxarg - std::log(std::exp(-maxarg) + std::exp(val - maxarg));
    }
    Real log_numer = LogSumExpOMP(logf_F, onesF);

    std::vector<Real> logf_B(sizeWB, 0.0);
    std::vector<Real> onesB(sizeWB, 1.0);
#pragma omp parallel for
    for (int i = 0; i < (int)sizeWB; i++) {
        Real val = -M + w_B[i] + DeltaF;
        Real maxarg = std::max(val, 0.0);

        logf_B[i] = -maxarg - std::log(std::exp(-maxarg) + std::exp(val - maxarg));
    }
    Real log_denom = LogSumExpOMP(logf_B, onesB);

    return log_numer - log_denom;
}

WhamTools::Real WhamTools::EXP(const std::vector<Real>& w_F)
{
    Real size = w_F.size();

    std::vector<Real> ones(size, 1.0);
    std::vector<Real> negw_F(size, 0.0);

    for (int i = 0; i < size; i++) {
        negw_F[i] = -w_F[i];
    }

    Real val = LogSumExpOMP(negw_F, ones);
    Real denom = std::log(size);

    return -(val - denom);
}

WhamTools::Real WhamTools::CalculateDeltaFBarIterative(const std::vector<Real>& w_F,
                                                       const std::vector<Real>& w_B,
                                                       int max_iterations, Real tol)
{
    Real DeltaF = 0.0;

    for (int i = 0; i < max_iterations; i++) {
        Real DeltaFold = DeltaF;
        DeltaF = DeltaFold - CalculateBAR(w_F, w_B, DeltaFold);

        Real scale = std::max(std::abs(DeltaFold), 1.0);
        Real relativeChange = std::abs(DeltaF - DeltaFold) / scale;

        if (relativeChange < tol) {
            break;
        }
    }

    return DeltaF;
}

WhamTools::Real WhamTools::CalculateDeltaFBarBisection(const std::vector<Real>& w_F,
                                                       const std::vector<Real>& w_B,
                                                       int max_iterations)
{
    ASSERT((max_iterations > 0), "The maximum number of iterations must be positive.");

    Real lowerB = -EXP(w_B);
    Real upperB = EXP(w_F);
    if (upperB < lowerB) {
        std::swap(upperB, lowerB);
    }

    Real FLowerB = CalculateBAR(w_F, w_B, lowerB);
    Real FUpperB = CalculateBAR(w_F, w_B, upperB);
    const Real tolerance = 1e-7;

    if (std::abs(FLowerB) < tolerance) {
        return lowerB;
    }
    if (std::abs(FUpperB) < tolerance) {
        return upperB;
    }

    for (int i = 0; i < max_iterations && FLowerB * FUpperB > 0.0; i++) {
        Real width = upperB - lowerB;
        if (width == 0.0) {
            width = 1.0;
        }

        lowerB -= width;
        upperB += width;
        FLowerB = CalculateBAR(w_F, w_B, lowerB);
        FUpperB = CalculateBAR(w_F, w_B, upperB);
    }

    ASSERT((FLowerB * FUpperB <= 0.0), "The initial guesses must bracket a root.");

    Real mid = 0.5 * (lowerB + upperB);
    for (int i = 0; i < max_iterations; i++) {
        mid = 0.5 * (lowerB + upperB);
        Real FMid = CalculateBAR(w_F, w_B, mid);

        Real scale = std::max(std::abs(mid), 1.0);
        if (std::abs(FMid) < tolerance || std::abs(upperB - lowerB) / scale < tolerance) {
            return mid;
        }

        if (FLowerB * FMid <= 0.0) {
            upperB = mid;
            FUpperB = FMid;
        } else {
            lowerB = mid;
            FLowerB = FMid;
        }
    }

    return mid;
}

WhamTools::Real WhamTools::Uwham_NLL_equation(const std::vector<Real>& f_k,
                                              const Matrix<Real>& BUki, const std::vector<Real>& N)
{
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    // Get the total N
    Real Ntot = VectorOP::VectorSum(N);

    // get the fraction of N/Ntot
    std::vector<Real> N_fraction(Nsim, 0);
    for (int i = 0; i < Nsim; i++) {
        N_fraction[i] = N[i] / Ntot;
    }

    ASSERT((f_k.size() == Nsim),
           "The dimension of fk does not match that of the number of simulation.");

    // Calculates the first part of the equation
    Real firstPart = 0.0;
    for (int i = 0; i < Nsim; i++) {
        firstPart += N[i] * f_k[i];
    }

    // Calculates the second part of the equation
    Real secondPart = 0.0;
#pragma omp parallel for reduction(+ : secondPart)
    for (int i = 0; i < Ndata; i++) {
        std::vector<Real> temp(Nsim);
        for (int j = 0; j < Nsim; j++) {
            temp[j] = f_k[j] - BUki(j, i);
        }

        secondPart += WhamTools::LogSumExp(temp, N_fraction);
    }

    return -firstPart + secondPart;
}
