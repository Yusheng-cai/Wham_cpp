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
    chosen_data_.resize(size_,dimension_);

    for (int i=0;i<size_;i++)
    {
        ASSERT((Totaldata_[i].size() >= larger_col_), "The inputted column is wrong, the total size of the column is " << Totaldata_[i].size() << \
        " while it is trying to access the " << larger_col_ << "th item.");

        int i_inner = i + skipFromBeginning_;
        for (int j=0;j<dimension_;j++)
        {
            chosen_data_(i,j) = Totaldata_[i_inner][columns_[j] - 1];
        }
    }

    for (int i=0;i<size_;i++)
    {
        for (int j=0;j<dimension_;j++)
        {
            std::cout << chosen_data_(i,j) << std::endl;
        }
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
            Mean_[j] += chosen_data_(i,j);
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