#include "Driver.h"

Driver::Driver(const ParameterPack& pack, const CommandLineArguments& cmd)
    : pack_(const_cast<ParameterPack&>(pack))
{
    // read in the absolute path of the timeseries input, if not read then timeseries will assume it
    // is the current dir
    std::string abspath;
    cmd.readString("abspath", CommandLineArguments::Keys::Optional, abspath);

    // find all the instances of the timeseries block
    auto timeSeriesPacks = pack_.findParamPacks("timeseries", ParameterPack::KeyType::Required);

    for (int i = 0; i < timeSeriesPacks.size(); i++) {
        TimeSeriesInputPack input = {const_cast<ParameterPack&>(*timeSeriesPacks[i]), abspath};
        timeSeries_.push_back(TimeSeriesPtr(new TimeSeries(input)));
    }

    initializeWham();
    initializeReweight();
    initializeTimeSeriesOperations();
}

void Driver::initializeTimeSeriesOperations()
{
    auto operationPacks = pack_.findParamPacks("tsoperation", ParameterPack::KeyType::Optional);

    if (operationPacks.size() != 0) {
        for (int i = 0; i < operationPacks.size(); i++) {
            std::string operationType;
            operationPacks[i]->ReadString("type", ParameterPack::KeyType::Required, operationType);

            TSInput input = {const_cast<ParameterPack&>(*operationPacks[i]), timeSeries_};

            timeSeriesOperations_.push_back(TimeSeriesOperationPtr(
                timeseriesOP::Factory::instance().create(operationType, input)));
        }
    }
}

void Driver::initializeWham()
{
    auto whamPack = pack_.findParamPacks("wham", ParameterPack::KeyType::Optional);

    if (whamPack.size() != 0) {
        for (int i = 0; i < whamPack.size(); i++) {
            std::string whamType;
            whamPack[i]->ReadString("type", ParameterPack::KeyType::Required, whamType);

            WhamInput input = {const_cast<ParameterPack&>(pack_), timeSeries_};
            whamCalculations_.push_back(
                WhamPtr(WhamRegistry::Factory::instance().create(whamType, input)));

            std::string whamName = whamCalculations_[i]->getName();
            auto it = whamNameToIndex_.find(whamName);

            ASSERT((it == whamNameToIndex_.end()),
                   "The name of wham " << whamName << " is registered twice.");
            whamNameToIndex_.insert(std::make_pair(whamName, i));
        }
    }
}

void Driver::initializeReweight()
{
    auto reweightPack = pack_.findParamPacks("Reweight", ParameterPack::KeyType::Optional);

    for (int i = 0; i < reweightPack.size(); i++) {
        std::string whamName;
        std::string reweighttype;

        // specify which WHAM this is acting on
        reweightPack[i]->ReadString("wham", ParameterPack::KeyType::Required, whamName);
        reweightPack[i]->ReadString("type", ParameterPack::KeyType::Required, reweighttype);

        // find the WHAM
        auto it = whamNameToIndex_.find(whamName);
        ASSERT((it != whamNameToIndex_.end()),
               "The name of wham " << whamName << " is not registered.");
        int index = it->second;

        ReweightInput input = {whamCalculations_[index].get(),
                               const_cast<ParameterPack&>(*reweightPack[i])};
        reweights_.push_back(
            ReweightPtr(ReweightRegistry::Factory::instance().create(reweighttype, input)));
    }
}

void Driver::calculate()
{
    // calculate timeseries first
    for (int i = 0; i < timeSeries_.size(); i++) {
        timeSeries_[i]->calculate();
    }

    for (int i = 0; i < whamCalculations_.size(); i++) {
        whamCalculations_[i]->calculate();
    }

    for (int i = 0; i < reweights_.size(); i++) {
        reweights_[i]->calculate();
    }

    for (int i = 0; i < timeSeriesOperations_.size(); i++) {
        timeSeriesOperations_[i]->calculate();
    }
}

void Driver::finishCalculate()
{
    for (int i = 0; i < whamCalculations_.size(); i++) {
        whamCalculations_[i]->finishCalculate();
    }
}

void Driver::printOutput()
{
    for (int i = 0; i < whamCalculations_.size(); i++) {
        whamCalculations_[i]->printOutput();
    }

    for (int i = 0; i < reweights_.size(); i++) {
        reweights_[i]->printOutput();
    }

    for (int i = 0; i < timeSeries_.size(); i++) {
        timeSeries_[i]->printOutput();
    }

    for (int i = 0; i < timeSeriesOperations_.size(); i++) {
        timeSeriesOperations_[i]->print();
    }
}
