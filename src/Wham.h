#pragma once
#include "Eigen/Dense"
#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"
#include "tools/Constants.h"
#include "tools/GenericFactory.h"
#include "Array.h"
#include "VectorOperations.h"
#include "Bias.h"
#include "Bin.h"

#include <vector>
#include <array>
#include <string>
#include <chrono>
#include <memory>
#include <functional>
#include <algorithm>
#include <map>
#include <numeric>

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
        using tsptr= std::shared_ptr<TimeSeries>;
        using Real = CommonTypes::Real;
        using valueFunction = std::function<void(std::string)>;

        Wham(const WhamInput& input);
        virtual ~Wham(){};

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

        std::string getName() {return name_;}
        int getDimension() const {return dimension_;}

        // Printing functions to be registered
        void printTimeSeriesBins(std::string name);

        // print out the force -dU/dx of the bias
        void printForce(std::string name);

        // print out the autocorrelation of the data 
        void printAutocorrelation(std::string name);

        // print out the average quantities of the data 
        void printAverage(std::string name);
    
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

        int precision_=3;

        // name of the wham  --> defaulted to "w"
        std::string name_ = "w";

        // histogram for each dimension of data
        std::vector<std::vector<std::vector<Real>>> histogram_;

        // The bins used in the calculation
        std::vector<Bin> Bins_;

        // The averages and standard deviations of the timeseries 
        std::vector<std::vector<Real>> Averages_;
        std::vector<std::vector<Real>> Std_;
};

namespace WhamRegistry
{
    using Base = Wham;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,const WhamInput&>;

    template<typename D>
    using registry = RegisterInFactory<Base, D, Key, const WhamInput&>;
};

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
