#pragma once

#include "Array.h"
#include "tools/CommonTypes.h"

#include <vector>

namespace WhamTools
{
    using Real = CommonTypes::Real;

    // Pass in vector is calculated as Log(\sum(N*exp(vector)))
    Real LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N);

    // Calculate log sum exp but with omp
    Real LogSumExpOMP(const std::vector<Real>& vector, const std::vector<Real>& N);

    // find the norm of a vector
    Real NormVector(const std::vector<Real>& vector);

    // find the hessian matrix of the UWham NLL equation
    Matrix<Real> Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    // find the gradient vector of the UWham NLL equation
    std::vector<Real> Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    // find the gradient vector of the Bwham NLL equation
    std::vector<Real> BGradient(const Matrix<Real>& BWil, const std::vector<Real>& Ml, const std::vector<Real>& N, \
    const std::vector<Real>& fk);

    // find the lnWi in UWham
    std::vector<Real> calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N);

    // find the lnpl in Bwham
    std::vector<Real> calculatelnpl(const Matrix<Real>& BWil, const std::vector<Real>& Ml, const std::vector<Real>& N, \
    const std::vector<Real>& fk);

    // Calculate the BAR value
    Real CalculateBAR(const std::vector<Real>& w_F, const std::vector<Real>& w_B, Real DeltaF);

    // Estimate free energy difference using one-sided (unidirectional) exponential averaging (EXP)
    Real EXP(const std::vector<Real>& w_F);

    // Estimate free energy difference using BAR --> using iterative method
    Real CalculateDeltaFBarIterative(const std::vector<Real>& w_F, const std::vector<Real>& w_B, int maxiterations=500, Real tolerance=1e-7);

    // Estimate free energy difference using BAR --> bisection method
    Real CalculateDeltaFBarBisection(const std::vector<Real>& w_F, const std::vector<Real>& w_B, int maxiterations=500);

    // Calculates the Uwham NLL equation
    Real Uwham_NLL_equation(const std::vector<Real>& f_k, const Matrix<Real>& BUji, const std::vector<Real>& N);
};
