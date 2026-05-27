#pragma once
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "Eigen/Dense"
#include "TimeSeries.h"
#include "WhamTools.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/Constants.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

struct WhamInput
{
    using tsptr = std::shared_ptr<TimeSeries>;
    ParameterPack& pack_;
    std::vector<tsptr>& VectorTimeSeries_;
};

class Wham
{
public:
    using Biasptr = std::unique_ptr<Bias>;
    using tsptr = std::shared_ptr<TimeSeries>;
    using Real = CommonTypes::Real;
    using valueFunction = std::function<void(std::string)>;

    Wham(const WhamInput& input);
    virtual ~Wham() {};

    void registerOutput(std::string name, valueFunction func);
    valueFunction& printOutputFromName(std::string name);

    virtual void calculate() = 0;
    virtual void initializeBias();
    virtual void initializeTimeSeries();
    void binTimeSeries();
    void initializeBins();

    // check if all the outputs are registered
    void isRegistered();

    virtual void printOutput();
    virtual void finishCalculate() {};

    virtual std::string type() = 0;

    std::string getName() { return name_; }
    int getDimension() const { return dimension_; }

    // Printing functions to be registered
    void printTimeSeriesBins(std::string name);

    // print out the force -dU/dx of the bias
    void printForce(std::string name);

    // print out the autocorrelation of the data
    void printAutocorrelation(std::string name);

    // print out the average quantities of the data
    void printAverage(std::string name);

    void printdataFE(std::string name);

protected:
    std::vector<tsptr>& VectorTimeSeries_;

    std::vector<Real> N_;

    // total data
    std::vector<std::vector<Real>> xi_;

    std::map<std::string, valueFunction> MapNameToFunction_;

    // output names as well as output file names
    std::vector<std::string> VectorOutputNames_;
    std::vector<std::string> VectorOutputFileNames_;

    // the parameter pack
    ParameterPack& pack_;
    ParameterPack* whamPack_;

    // the vector of all the bias
    std::vector<Biasptr> Biases_;

    // The dimension of the timeseries
    std::vector<int> dimensions_;
    int dimension_;

    // Total number of data
    int Ntot_ = 0;

    int precision_ = 3;

    bool verbose_ = false;

    // name of the wham  --> defaulted to "w"
    std::string name_ = "w";

    // histogram for each dimension of data
    std::vector<std::vector<std::vector<Real>>> histogram_;

    // normalized histogram --> -ln(pk)
    std::vector<std::map<std::vector<int>, Real>> dataFE_;

    // The bins used in the calculation
    std::vector<Bin> Bins_;

    // The averages and standard deviations of the timeseries
    std::vector<std::vector<Real>> Averages_;
    std::vector<std::vector<Real>> Std_;

    // Combined input
    bool combined_input_ = false;
};

namespace WhamRegistry {
using Base = Wham;
using Key = std::string;

using Factory = GenericFactory<Base, Key, const WhamInput&>;

template <typename D> using registry = RegisterInFactory<Base, D, Key, const WhamInput&>;
}; // namespace WhamRegistry
