#pragma once
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/FileSystem.h"
#include "DataFileParser.h"
#include "Array.h"
#include "FFT.h"

#include <vector>
#include <array>
#include <map>
#include <string>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cmath>
#include <iostream>
#include <functional>


struct TimeSeriesInputPack
{
    ParameterPack& pack_;
    std::string abspath_;
};

class TimeSeries
{
    public:
        using Real = CommonTypes::Real;
        using ComplexReal = std::complex<Real>;
        using Iterator = std::vector<std::vector<Real>>::iterator;
        using cIterator = std::vector<std::vector<Real>>::const_iterator;
        using valueFunction = std::function<void(std::string)>;
        using calcFunction = std::function<void()>;

        TimeSeries(const TimeSeriesInputPack& input);
        ~TimeSeries(){};


        // getters
        int getDimension() const {return dimension_;}

        // size of the data (number of data)
        int getSize() const {return size_;}

        void calculate();
        void calculateAutoCorrelation();
 
        // find the mean and variance of the TimeSeries
        void findMean();
        void findVar();

        // find normalized data once mean and variance has been found  
        void findNormalizedData();

        // print output
        void printOutput();

        // register outputs to print
        void registerOutputFunctions(std::string name, valueFunction printFunction);

        // register calculation to print
        void registerCalculateFunctions(std::string name, calcFunction calFunc);

        // check if output is Valid
        void checkOutputValidity();

        // read the chosen data based on user input from the total data 
        void readChosenData();

        void printAC(std::string name);
        valueFunction& printOutputFromName(std::string outputName);
        calcFunction& calculateFromName(std::string calcName);

        // get the data raw pointer underneath
        std::vector<Real>* data() {return chosen_data_.data();}
        Iterator begin() {return chosen_data_.begin();}
        Iterator end() {return chosen_data_.end();}
        cIterator cbegin() {return chosen_data_.cbegin();}
        cIterator cend() {return chosen_data_.cend();}
        std::vector<Real>& operator[](int i) { return chosen_data_[i];}

        // get the mean of the timeseries 
        std::vector<Real> getMean() const {return Mean_;}
        std::vector<Real> getstd() const {return std_;}

        // get the lag time / autocorrelation time 
        std::vector<Real> getLagTime() const {return lag_time_;}

    
    private:
        std::string path_;
        DataFileParser parser;

        // entire data
        //Matrix<Real> chosen_data_;
        std::vector<std::vector<Real>> chosen_data_;
        std::vector<std::vector<Real>> Totaldata_; 

        // Sometimes it is ideal to skip some data in the beginning because the simulation has not reached eq. yet
        int skipFromBeginning_ = 0;

        // vector of columns
        std::vector<int> columns_;

        // size of the data
        int size_;

        // dimension of the data
        int dimension_;

        // the larger column of the inputs, e.g [ 1, 2, 3 ] -> 3
        int larger_col_;

        // the number of data to skip 
        int skipevery_=0;

        // Find mean for each dimension of the timeseries
        std::vector<Real> Mean_, Variance_, std_;


        // vector of autocorrelations
        std::vector<std::vector<Real>> AC_vector_;
        std::vector<Real> lag_time_;

        // normalized data 
        std::vector<std::vector<Real>> normalized_Data_;

        std::vector<std::string> outputNames_;
        std::vector<std::string> outputFileNames_;

        // Map from name to function
        std::map<std::string, valueFunction> MapNameToOutput_;

        // Map from name to calculate function
        std::map<std::string, calcFunction> MapNameToCalculate_;
};