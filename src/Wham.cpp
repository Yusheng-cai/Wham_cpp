#include "Wham.h"

Wham::Wham(const ParameterPack& pack)
{
    // find all the instances of the timeseries block
    auto TsPacks = pack.findParamPacks("timeseries", ParameterPack::KeyType::Required);

    for (int i= 0 ;i < TsPacks.size();i++)
    {
        VectorTimeSeries_.push_back(TimeSeries(*TsPacks[i]));
    }
}

WhamTools::Real WhamTools::LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N)
{
    ASSERT((vector.size() == N.size()), "The size of the vector is not equal to the size of N.");

    // Find the max of the vector
    auto it = std::max_element(vector.begin(), vector.end());
    Real maxVal = *it;

    Real sum = 0.0;
    for (int i=0;i<vector.size();i++)
    {
        sum += N[i] * std::exp(vector[i]-maxVal);
    }

    sum = std::log(sum) * maxVal;

    return sum;
}