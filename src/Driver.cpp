#include "Driver.h"

Driver::Driver(const ParameterPack& pack, const CommandLineArguments& cmd)
{
    // read in the absolute path of the timeseries input, if not read then timeseries will assume it is the current dir
    std::string abspath;
    bool read = cmd.readString("abspath", CommandLineArguments::Keys::Optional,abspath);

    // find all the instances of the timeseries block
    auto TsPacks = pack.findParamPacks("timeseries", ParameterPack::KeyType::Required);

    for (int i= 0 ;i < TsPacks.size();i++)
    {
        TimeSeriesInputPack input = { const_cast<ParameterPack&>(*TsPacks[i]), abspath};
        VectorTimeSeries_.push_back(tsptr(new TimeSeries(input)));
    }

    InitializeWham(pack);
    InitializeTSoperation(pack);
}

void Driver::InitializeTSoperation(const ParameterPack& pack)
{
    auto TSPack = pack.findParamPacks("tsoperation", ParameterPack::KeyType::Optional);

    if (TSPack.size() != 0)
    {
        for (int i=0;i<TSPack.size();i++)
        {
            std::string type_;
            TSPack[i] -> ReadString("type", ParameterPack::KeyType::Required,type_);

            TSInput input = {const_cast<ParameterPack&>(*TSPack[i]), VectorTimeSeries_};

            VectorTimeSeriesOP_.push_back(TSopptr(timeseriesOP::Factory::instance().create(type_, input)));
        }
    }
}

void Driver::InitializeWham(const ParameterPack& pack)
{
    auto whamPack = pack.findParamPacks("wham", ParameterPack::KeyType::Optional);

    if (whamPack.size() != 0)
    {
        for (int i=0;i<whamPack.size();i++)
        {
            std::string whamType;
            whamPack[i]->ReadString("type", ParameterPack::KeyType::Required,whamType);

            WhamInput input = {const_cast<ParameterPack&>(pack), VectorTimeSeries_};
            VectorWhamCalc_.push_back(Whamptr(WhamRegistry::Factory::instance().create(whamType, input)));
        }
    }
}
    
void Driver::calculate()
{
    // calculate timeseries first 
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        VectorTimeSeries_[i] -> calculate();
    }

    for (int i=0;i<VectorWhamCalc_.size();i++)
    {
        VectorWhamCalc_[i] -> calculate();
    }

    for (int i=0;i<VectorTimeSeriesOP_.size();i++)
    {
        VectorTimeSeriesOP_[i] -> calculate();
    }

}

void Driver::finishCalculate()
{
    for (int i=0;i<VectorWhamCalc_.size();i++)
    {
        VectorWhamCalc_[i] -> finishCalculate();
    }
}

void Driver::printOutput()
{
    for (int i=0;i<VectorWhamCalc_.size();i++)
    {
        VectorWhamCalc_[i] -> printOutput();
    }

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        VectorTimeSeries_[i] -> printOutput();
    }

    for (int i=0;i<VectorTimeSeriesOP_.size();i++)
    {
        VectorTimeSeriesOP_[i] -> print();
    }
}