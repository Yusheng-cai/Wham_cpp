#include "TSoperation.h"

TSoperation::TSoperation(const TSInput& input)
:VectorTS_(input.vectorTS_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, VectorOutputFileNames_);

    ASSERT((VectorOutputNames_.size() == VectorOutputFileNames_.size()), "Output files does not match output names.");

    outputs_ = OutputFuncPtr(new Output());

    combineData();

    outputs_ -> registerOutputFunc("totaldata", [this](std::string name) -> void {printTotalData(name);});
    outputs_ -> registerOutputFunc("totaldatalength", [this](std::string name) -> void {printTotalDataLength(name);});
    outputs_ -> registerOutputFunc("averages", [this](std::string name) -> void {printMean(name);});
    outputs_ -> registerOutputFunc("lagtime", [this](std::string name) -> void {printLagTime(name);});
}

void TSoperation::print()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        outputs_ -> getOutputFuncByName(VectorOutputNames_[i])(VectorOutputFileNames_[i]);
    }
}

void TSoperation::printLagTime(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# TimeSeries LagTime(timestep)\n";

    for (int i=0;i<LongestLagTimes_.size();i++)
    {
        ofs << i + 1 << " " << LongestLagTimes_[i] << "\n";
    }

    ofs.close();
}

void TSoperation::printMean(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int index = 1;
    ofs << "# ";

    for (int i=0;i<averages_[0].size();i++)
    {
        ofs << "OP" << i+1 << " ";
    }

    for (int i=0;i<averages_[0].size();i++)
    {
        ofs << "std" << i+1 << " ";
    }
    ofs << "\n";

    for (int i=0;i<averages_.size();i++)
    {
        ofs << i+1 << " ";
        for (auto num : averages_[i])
        {
            ofs << num << " ";
        }

        for (auto num : std_[i])
        {
            ofs << num << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void TSoperation::printTotalDataLength(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<TotalDataLength_.size();i++)
    {
        ofs << TotalDataLength_[i] << " ";
    }

    ofs.close();
}

void TSoperation::combineData()
{
    for (int i=0;i<VectorTS_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTS_[i]->begin(), VectorTS_[i]->end());
        TotalDataLength_.push_back(VectorTS_[i]->getSize());

        averages_.push_back(VectorTS_[i]->getMean());
        std_.push_back(VectorTS_[i]->getstd());
        LongestLagTimes_.push_back(VectorTS_[i] -> getLongestLagTime());
    }
}

void TSoperation::printTotalData(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<xi_.size();i++)
    {
        for (int j=0;j<xi_[i].size();j++)
        {
            ofs << xi_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}