#include "TSoperation.h"

TSoperation::TSoperation(const TSInput& input)
:VectorTS_(input.vectorTS_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    pack_.ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);

    outputs_ = OutputFuncPtr(new Output());

    combineData();

    outputs_ -> registerOutputFunc("totaldata", [this](std::string name) -> void {printTotalData(name);});
    outputs_ -> registerOutputFunc("totaldatalength", [this](std::string name) -> void {printTotalDataLength(name);});
}

void TSoperation::print()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        outputs_ -> getOutputFuncByName(VectorOutputNames_[i])(VectorOutputFileNames_[i]);
    }
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
    }
}

void TSoperation::printTotalData(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    std::cout << "Xi size = " << xi_.size() << std::endl;

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