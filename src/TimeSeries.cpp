#include "TimeSeries.h"

TimeSeries::TimeSeries(const ParameterPack& pack)
{
    pack.ReadString("path", ParameterPack::KeyType::Required, path_);
    pack.ReadNumber("skipfrombeginning", ParameterPack::KeyType::Optional, skipFromBeginning_);
    pack.ReadVectorNumber("columns", ParameterPack::KeyType::Required, columns_);

    auto it = std::max_element(columns_.begin(), columns_.end());
    larger_col_ = *it;

    // find out the dimension of the time series
    dimension_ = columns_.size();

    // Have the parser parse the inputted file
    parser.ParseFile(path_, Totaldata_);

    // Find out the total size of the data
    size_ = Totaldata_.size() - skipFromBeginning_;

    // Resize the chosen data accordingly
    chosen_data_.resize(size_);
    int index=0;
    for (auto it = Totaldata_.begin()+skipFromBeginning_;it != Totaldata_.end();it++)
    {
        std::vector<Real> temp;
        temp.resize(dimension_);
        for (int i=0;i<dimension_;i++)
        {
            temp[i] = (*it)[columns_[i]-1];
        }
        chosen_data_[index] = temp;
        index ++;
    }

    findMean(); 
}

void TimeSeries::findMean()
{
    Mean_.resize(dimension_);

    for (int i=0;i<dimension_;i++)
    {
        Mean_[i] = 0.0;
    }

    for (int i=0;i<size_;i++)
    {
        for(int j=0;j<dimension_;j++)
        {
            Mean_[j] += chosen_data_[i][j];
        }
    }

    // Find the mean of the system
    for (int i=0;i<dimension_;i++)
    {
        Mean_[i] = Mean_[i]/size_;
    }
}

void TimeSeries::findVar()
{

}