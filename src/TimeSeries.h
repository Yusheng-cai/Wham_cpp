#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "DataFileParser.h"
#include "Array.h"

#include <vector>
#include <array>
#include <string>
#include <algorithm>


class TimeSeries
{
    public:
        using Real = CommonTypes::Real;

        TimeSeries(const ParameterPack& pack);
        ~TimeSeries(){};

        // getters
        int getDimension() const {return dimension_;}
        int getSize() const {return size_;}

        
        // find the mean and variance of the TimeSeries
        void findMean();
        void findVar();

    
    private:
        std::string path_;
        DataFileParser parser;

        // entire data
        Matrix<Real> chosen_data_;
        std::vector<std::vector<Real>> Totaldata_; 

        // Sometimes it is ideal to skip some data
        int skip_ = 1;

        // vector of columns
        std::vector<int> columns_;

        // size of the data
        int size_;

        // dimension of the data
        int dimension_;

        // the larger column of the inputs, e.g [ 1, 2, 3 ] -> 3
        int larger_col_;

        // Find mean for each dimension of the timeseries
        std::vector<Real> Mean_, Variance_;
};