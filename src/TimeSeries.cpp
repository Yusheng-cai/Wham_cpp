#include "TimeSeries.h"

TimeSeries::TimeSeries(const TimeSeriesInputPack& input)
{
    input.pack_.ReadString("path", ParameterPack::KeyType::Required, path_);
    input.pack_.ReadNumber("skipfrombeginning", ParameterPack::KeyType::Optional, skipFromBeginning_);
    input.pack_.ReadVectorNumber("columns", ParameterPack::KeyType::Required, columns_);
    input.pack_.ReadVectorNumber("AC_dimension", ParameterPack::KeyType::Optional,AC_dimensions_);


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
    std::cout << "Reading data from file " <<path_ << " data size = " << Totaldata_.size() << std::endl;

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

void TimeSeries::calculate()
{
    if (AC_dimensions_.size() > 0)
    {
        for (int i=0;i<AC_dimensions_.size();i++)
        {
            int dim = AC_dimensions_[i] - 1;

            // to calculate autocorrelation, we must split the data up
            std::vector<Real> data_(chosen_data_.size()*2,0.0);
            int datasize = chosen_data_.size();
            int hdatasize = datasize/2;
            int otherhalf = datasize - hdatasize;

            for (int j=0;j<hdatasize;j++)
            {
                int indexCS = datasize - hdatasize + j - 1;
                data_[j] = chosen_data_[indexCS][dim];
            }

            for (int j=0;j<otherhalf;j++)
            {
                int indexd = hdatasize + datasize + j - 1;
                data_[indexd] = chosen_data_[j][dim]; 
            }

            std::vector<Real> AC;
            calculateAutoCorrelation(data_,AC);

            for (int i=0;i<AC.size();i++)
            {
                std::cout << AC[i] << std::endl;
            }
        }
    }
}

void TimeSeries::calculateAutoCorrelation(const std::vector<Real>& data, std::vector<Real>& AC)
{
    int datasize = data.size();

    std::vector<ComplexReal> fftComplex_;
    std::vector<ComplexReal> input_(datasize);

    for (int i=0;i<datasize;i++)
    {
        ComplexReal number(data[i],0);
        input_[i] = number;
    }

    FFT::fft(input_, input_);

    std::vector<ComplexReal> squared(datasize);
    std::cout << "Made squared" << std::endl;
    for (int i=0;i<datasize;i++)
    {
        squared[i] = std::pow(fftComplex_[i].real(),2.0) + std::pow(fftComplex_[i].imag(),2.0);
    }
    std::cout << "Done squared" << std::endl;

    std::vector<ComplexReal> output_;
    FFT::fft(squared, output_);

    for (int i=0;i<datasize;i++)
    {
        AC[i] = output_[i].real()/datasize;
    }
}

