#include "TimeSeries.h"

TimeSeries::TimeSeries(const TimeSeriesInputPack& input)
{
    input.pack_.ReadString("path", ParameterPack::KeyType::Required, path_);
    input.pack_.ReadNumber("skipfrombeginning", ParameterPack::KeyType::Optional, skipFromBeginning_);
    input.pack_.ReadVectorNumber("columns", ParameterPack::KeyType::Required, columns_);
    input.pack_.ReadVectorNumber("AC_dimension", ParameterPack::KeyType::Optional,AC_dimensions_);
    bool output = input.pack_.ReadString("output", ParameterPack::KeyType::Optional, OutputName_);

    if (output)
    {
        ofs_.open(OutputName_);
    }


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

    // find the mean of the data set
    findMean();

    // find the variance of the data size
    findVar();

    // normalize the data 
    findNormalizedData();
}

void TimeSeries::findNormalizedData()
{
    std::vector<Real> zeros(dimension_,0.0);

    int Numdata = chosen_data_.size();
    normalized_Data_.resize(Numdata, zeros);

    for (int i=0;i<Numdata;i++)
    {
        for (int j=0;j<dimension_;j++)
        {
            normalized_Data_[i][j] = (chosen_data_[i][j] - Mean_[j])/std_[j];
        }
    }
}

void TimeSeries::findMean()
{
    Mean_.resize(dimension_,0.0);

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
    // find the number of data 
    int NumData = chosen_data_.size();

    // resize variance 
    Variance_.resize(dimension_,0.0);
    std_.resize(dimension_,0.0);

    for (int i=0;i<NumData;i++)
    {
        for (int j=0;j<dimension_;j++)
        {
            Real diff = chosen_data_[i][j] - Mean_[j];
            Variance_[j] += std::pow(diff,2);
        }
    }

    for (int j=0;j<dimension_;j++)
    {
        Variance_[j] /= NumData;
        std_[j] = std::sqrt(Variance_[j]);
    }
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
            
            for (int j=0;j<otherhalf;j++)
            {
                int indexCS = hdatasize + j;
                data_[j] = normalized_Data_[indexCS][dim];
            }

            for (int j=0;j<hdatasize;j++)
            {
                int indexd = hdatasize + datasize + j + 1;
                data_[indexd] = normalized_Data_[j][dim]; 
            }


            calculateAutoCorrelation(data_,AC_);
        }
    }
}

void TimeSeries::calculateAutoCorrelation(const std::vector<Real>& data, std::vector<Real>& AC)
{
    int datasize = data.size();
    AC.clear();
    AC.resize(datasize,0.0);

    std::vector<ComplexReal> fftComplex_;
    std::vector<ComplexReal> input_(datasize);

    for (int i=0;i<datasize;i++)
    {
        ComplexReal number(data[i],0);
        input_[i] = number;
    }

    FFT::fft(input_, fftComplex_);

    std::vector<ComplexReal> squared(datasize);
    for (int i=0;i<datasize;i++)
    {
        Real square = std::pow(fftComplex_[i].real(),2.0) + std::pow(fftComplex_[i].imag(),2.0);
        ComplexReal number(square,0.0);
        squared[i] = number;
    }

    std::vector<ComplexReal> output_;
    FFT::ifft(squared, output_);

    for (int i=0;i<datasize;i++)
    {
        AC[i] = 2*output_[i].real()/datasize;
    }
}

void TimeSeries::printOutput()
{
    if (ofs_.is_open())
    {
        int datasize = chosen_data_.size();

        for (int i=0;i<datasize;i++)
        {
            ofs_ << AC_[i] << "\n";
        }

        ofs_.close();
    }
}