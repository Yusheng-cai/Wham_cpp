#pragma once
#include "Assert.h"
#include "CommonTypes.h"

#include <array>
#include <functional>
#include <map>
#include <vector>

class Output
{
public:
    using outputFunc = std::function<void(std::string name)>;

    Output() = default;

    void registerOutputFunc(std::string name, outputFunc func);
    outputFunc& getOutputFuncByName(std::string name);

    bool isregistered(std::string name);

private:
    std::map<std::string, outputFunc> MapNameToOutputFunc_;
};