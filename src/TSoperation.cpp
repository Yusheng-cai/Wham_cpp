#include "TSoperation.h"

TSoperation::TSoperation(const TSInput& input)
:VectorTS_(input.vectorTS_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    pack_.ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);

    outputs_ = OutputFuncPtr(new Output());

    combineData();

    outputs_ -> registerOutputFunc("totaldata", [this](std::string name) -> void {printTotalData(name);});
}

void TSoperation::print()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        outputs_ -> getOutputFuncByName(VectorOutputNames_[i])(VectorOutputFileNames_[i]);
    }
}

void TSoperation::combineData()
{
    for (int i=0;i<VectorTS_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTS_[i]->begin(), VectorTS_[i]->end());
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