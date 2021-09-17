#pragma once
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/FileSystem.h"
#include "DataFileParser.h"
#include "Array.h"
#include "FFT.h"

#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cmath>
#include <iostream>


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

        TimeSeries(const TimeSeriesInputPack& input);
        ~TimeSeries(){};


        // getters
        int getDimension() const {return dimension_;}
        int getSize() const {return size_;}

        void calculate();
        void calculateAutoCorrelation(const std::vector<Real>& data, std::vector<Real>& AC);
 
        // find the mean and variance of the TimeSeries
        void findMean();
        void findVar();

        // find normalized data once mean and variance has been found  
        void findNormalizedData();

        // print output
        void printOutput();

        // get the data raw pointer underneath
        std::vector<Real>* data() {return chosen_data_.data();}
        Iterator begin() {return chosen_data_.begin();}
        Iterator end() {return chosen_data_.end();}
        cIterator cbegin() {return chosen_data_.cbegin();}
        cIterator cend() {return chosen_data_.cend();}

    
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

        // Find mean for each dimension of the timeseries
        std::vector<Real> Mean_, Variance_, std_;

        // read autocorrelation dimensions
        std::vector<int> AC_dimensions_;

        // normalized data 
        std::vector<std::vector<Real>> normalized_Data_;

        std::string OutputName_;
        std::ofstream ofs_;

        std::vector<Real> AC_;
};