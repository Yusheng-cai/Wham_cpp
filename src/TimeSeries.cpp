#include "TimeSeries.h"

TimeSeries::TimeSeries(const TimeSeriesInputPack& input)
{
    // register output functions
    registerOutputFunctions("autocorrelation", [this](std::string name) -> void {printAC(name);});

    // register calculate functions
    registerCalculateFunctions("autocorrelation", [this]() -> void {calculateAutoCorrelation();});

    input.pack_.ReadString("path", ParameterPack::KeyType::Required, path_);
    input.pack_.ReadNumber("skipfrombeginning", ParameterPack::KeyType::Optional, skipFromBeginning_);
    input.pack_.ReadVectorNumber("columns", ParameterPack::KeyType::Required, columns_);
    input.pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional,outputNames_);
    input.pack_.ReadNumber("skip", ParameterPack::KeyType::Optional, skipevery_);
    skipevery_++;
    checkOutputValidity();
    input.pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, outputFileNames_);

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

    readChosenData();

    // find the mean of the data set
    findMean();

    // find the variance of the data size
    findVar();

    // normalize the data 
    findNormalizedData();

    // calculate the autocorrelation
    calculateAutoCorrelation();
}

void TimeSeries::readChosenData()
{
    // Have the parser parse the inputted file
    parser.ParseFile(path_, Totaldata_);
    
    // Find out the total size of the data
    std::cout << "Reading data from file " <<path_ << " data size = " << Totaldata_.size() << std::endl;

    // Resize the chosen data accordingly
    int index=0;

    if (skipFromBeginning_ < 0)
    {
        ASSERT(((skipFromBeginning_ + Totaldata_.size()) >= 0), "Total data size is " << Totaldata_.size() << " while you specified skipfrombeginning = " << skipFromBeginning_);
        skipFromBeginning_ = Totaldata_.size() + skipFromBeginning_;
    }

    int originalIndex = skipFromBeginning_;


    while (originalIndex < Totaldata_.size())
    {
        std::vector<Real> temp;
        temp.resize(dimension_);
        for (int j=0;j<dimension_;j++)
        {
            temp[j] = Totaldata_[originalIndex][columns_[j]-1];
        }
        chosen_data_.push_back(temp);
        index ++;
        originalIndex += skipevery_;
    }

    size_ = chosen_data_.size();
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
    for (int i=0;i<outputNames_.size();i++)
    {
        calculateFromName(outputNames_[i])();
    }
}

void TimeSeries::calculateAutoCorrelation()
{
    // we have an autocorrelation for each dimension
    AC_vector_.resize(dimension_);

    int N = chosen_data_.size();
    int N2 = N*2;

    for (int i=0;i<dimension_;i++)
    {
        // to calculate autocorrelation, we must split the data up
        std::vector<Real> data(N2,0.0);

        int hdatasize = N/2;
        int otherhalf = N - hdatasize;

        // resize the ith dimension AC to be data size large
        AC_vector_[i].resize(N);
        
        for (int j=0;j<otherhalf;j++)
        {
            int indexCS = hdatasize + j;
            data[j] = normalized_Data_[indexCS][i];
        }

        for (int j=0;j<hdatasize;j++)
        {
            int indexd = hdatasize + N + j + 1;
            data[indexd] = normalized_Data_[j][i]; 
        }

        std::vector<ComplexReal> fft;
        std::vector<ComplexReal> input(N2);

        for (int j=0;j<N2;j++)
        {
            ComplexReal number(data[j],0);
            input[j] = number;
        }

        FFT::fft(input, fft);

        std::vector<ComplexReal> squared(N2);
        for (int j=0;j<N2;j++)
        {
            Real square = std::pow(fft[j].real(),2.0) + std::pow(fft[j].imag(),2.0);
            ComplexReal number(square,0.0);
            squared[j] = number;
        }

        std::vector<ComplexReal> ifft;
        FFT::ifft(squared, ifft);

        for (int j=0;j<N;j++)
        {
            AC_vector_[i][j] = ifft[j].real()/N;
        }
    }

    lag_time_.resize(dimension_,0.0);
    for (int i=0;i<dimension_;i++)
    {
        for (int j=0;j<AC_vector_[i].size();j++)
        {
            if (AC_vector_[i][j] < 0.0)
            {
                break;
            }
            else
            {
                lag_time_[i] += AC_vector_[i][j];
            }
        }
    }

    for (int i=0;i<dimension_;i++)
    {
        lag_time_[i] = 1 + 2 * lag_time_[i];
    }

    if (dimension_ > 1)
    {
        auto long_it = std::max_element(lag_time_.begin(), lag_time_.end());
        longest_lag_time_ = *long_it;
    }
    else
    {
        longest_lag_time_ = lag_time_[0];
    }

    // get number of independent points
    numIndependentPoints_ = (int)(chosen_data_.size() / longest_lag_time_);
}

std::vector<std::vector<TimeSeries::Real>> TimeSeries::getIndependentsample()
{
    std::vector<int> Index = RandomTools::RandomPermute(chosen_data_.size(), numIndependentPoints_);

    std::vector<std::vector<Real>> data_independent(Index.size());

    for (int i=0;i<Index.size();i++)
    {
        data_independent[i] = chosen_data_[Index[i]];
    }

    return data_independent;
}

void TimeSeries::checkOutputValidity()
{
    if (outputNames_.size() != 0)
    {
        for (int i=0;i<outputNames_.size();i++)
        {
            std::string name = outputNames_[i];
            auto it = MapNameToOutput_.find(name);

            ASSERT((it != MapNameToOutput_.end()), "The output with name " << name << " is not found within context of timeseries.");
        }
    }
}

void TimeSeries::printOutput()
{
    for (int i=0;i<outputNames_.size();i++)
    {
        std::string name = outputNames_[i];

        printOutputFromName(name)(outputFileNames_[i]);
    }
}

void TimeSeries::registerOutputFunctions(std::string name, valueFunction function)
{
    auto it = MapNameToOutput_.find(name);

    ASSERT((it == MapNameToOutput_.end()), "The output " << name << " already in within Timeseries output.");

    MapNameToOutput_.insert(std::make_pair(name, function));
}

void TimeSeries::registerCalculateFunctions(std::string name, calcFunction calfunc)
{
    auto it = MapNameToCalculate_.find(name);
    ASSERT((it == MapNameToCalculate_.end()), "The calculation " << name << " already exist within Timeseries.");

    MapNameToCalculate_.insert(std::make_pair(name, calfunc));
}

TimeSeries::valueFunction& TimeSeries::printOutputFromName(std::string outputName)
{
    auto it = MapNameToOutput_.find(outputName);

    ASSERT((it != MapNameToOutput_.end()), "The output " << outputName << " is not registered.");

    return it -> second;
}

TimeSeries::calcFunction& TimeSeries::calculateFromName(std::string calcName)
{
    auto it = MapNameToCalculate_.find(calcName);

    ASSERT((it != MapNameToCalculate_.end()), "The calculation " << calcName << " is not registered.");

    return it -> second;
}

void TimeSeries::printAC(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);


    int numData = chosen_data_.size();
    ofs << "#";
    for (int i=0;i<dimension_;i++)
    {
        ofs <<  i+1 << " ";
    }
    ofs << "\n";

    for (int i=0;i<numData;i++)
    {
        for (int j=0;j<dimension_;j++)
        {
            ofs << AC_vector_[j][i] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}
