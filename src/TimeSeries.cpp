#include "TimeSeries.h"

TimeSeries::TimeSeries(const TimeSeriesInputPack& input)
{
    input.pack_.ReadString("path", ParameterPack::KeyType::Required, path_);
    input.pack_.ReadNumber("skipfrombeginning", ParameterPack::KeyType::Optional, skipFromBeginning_);
    input.pack_.ReadVectorNumber("columns", ParameterPack::KeyType::Required, columns_);


    if (input.abspath_.empty())
    {
        path_ = FileSystem::joinPath(FileSystem::getCurrentPath(), path_); 
    }
    else
    {
        path_ = FileSystem::joinPath(input.abspath_, path_);
    }

    for (int i=0;i<columns_.size();i++)
    {
        ASSERT((columns_[i] > 0), "The column is 1-based counting, it cannot be less than or equal to 0.");
    }

    auto it = std::max_element(columns_.begin(), columns_.end());
    larger_col_ = *it;

    // find out the dimension of the time series
    dimension_ = columns_.size();

    // Have the parser parse the inputted file
    parser.ParseFile(path_, Totaldata_);
    
    // Find out the total size of the data
    size_ = Totaldata_.size() - skipFromBeginning_;
    std::cout << "Total data size = " << Totaldata_.size() << std::endl;

    // Resize the chosen data accordingly
    chosen_data_.resize(size_);
    int index=0;
    for (int i = skipFromBeginning_;i < Totaldata_.size();i++)
    {
        std::vector<Real> temp;
        temp.resize(dimension_);
        for (int j=0;j<dimension_;j++)
        {
            temp[j] = Totaldata_[i][columns_[j]-1];
        }
        chosen_data_[index] = temp;
        index ++;
    }

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