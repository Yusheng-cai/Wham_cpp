#pragma once

#include <vector>
#include <numeric>
#include <algorithm>

namespace VectorOP
{
    template<typename T>
    T VectorSum(const std::vector<T>& vec);

    template<typename T>
    T Min(const std::vector<T>& vec);

    template<typename T>
    T Max(const std::vector<T>& vec);
};

template<typename T>
T VectorOP::VectorSum(const std::vector<T>& vec)
{
    T val = std::accumulate(vec.begin(), vec.end(),0);

    return val;
}

template<typename T>
T VectorOP::Min(const std::vector<T>& vec)
{
    auto minIt = std::min_element(vec.begin(), vec.end());
    return *minIt;
}

template<typename T>
T VectorOP::Max(const std::vector<T>& vec)
{
    auto maxIt = std::max_element(vec.begin(), vec.end());

    return *maxIt;
}