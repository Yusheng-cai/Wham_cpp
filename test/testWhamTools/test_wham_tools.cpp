#include "src/WhamTools.h"
#include "src/BwhamLBFGS.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
bool near(double actual, double expected, double tolerance)
{
    return std::abs(actual - expected) <= tolerance;
}

int require_near(const std::string& name, double actual, double expected, double tolerance)
{
    if (!near(actual, expected, tolerance))
    {
        std::cerr << name << " expected " << expected << " but got " << actual << "\n";
        return 1;
    }

    return 0;
}

int require_finite(const std::string& name, double actual)
{
    if (!std::isfinite(actual))
    {
        std::cerr << name << " expected a finite value but got " << actual << "\n";
        return 1;
    }

    return 0;
}
}

int main()
{
    int failures = 0;

    failures += require_near(
        "LogSumExp",
        WhamTools::LogSumExp(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}),
        std::log(2.0),
        1e-12);

    failures += require_near(
        "EXP",
        WhamTools::EXP(std::vector<double>{0.0, 0.0}),
        0.0,
        1e-12);

    failures += require_near(
        "CalculateDeltaFBarIterative",
        WhamTools::CalculateDeltaFBarIterative(
            std::vector<double>{0.2, 0.4, 0.6},
            std::vector<double>{-0.1, 0.1, 0.3},
            500,
            1e-12),
        0.15,
        1e-8);

    failures += require_near(
        "CalculateDeltaFBarBisection",
        WhamTools::CalculateDeltaFBarBisection(
            std::vector<double>{0.0, 0.0},
            std::vector<double>{0.0, 0.0},
            500),
        0.0,
        1e-8);

    failures += require_near(
        "CalculateDeltaFBarBisection nonzero",
        WhamTools::CalculateDeltaFBarBisection(
            std::vector<double>{0.2, 0.4, 0.6},
            std::vector<double>{-0.1, 0.1, 0.3},
            500),
        0.15,
        1e-6);

    Matrix<double> BWil(2, 3);
    for (int row=0; row<BWil.getNR(); row++)
    {
        for (int col=0; col<BWil.getNC(); col++)
        {
            BWil(row, col) = 0.0;
        }
    }
    std::vector<double> N{2.0, 2.0};
    std::vector<double> Ml{1.0, 0.0, 1.0};
    BwhamNLLInput input{BWil, N, Ml};
    BwhamNLL nll(input);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(2);
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(2);
    double value = nll(x, grad);

    failures += require_finite("BwhamNLL zero-count bin objective", value);
    for (int i=0; i<grad.size(); i++)
    {
        failures += require_finite("BwhamNLL zero-count bin gradient", grad[i]);
    }

    return failures == 0 ? 0 : 1;
}
