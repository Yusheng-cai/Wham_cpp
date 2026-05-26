#include "src/Wham.h"

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

    return failures == 0 ? 0 : 1;
}
