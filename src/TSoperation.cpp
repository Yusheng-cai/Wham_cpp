#include "TSoperation.h"

TSoperation::TSoperation(const TSInput& input)
:VectorTS_(input.vectorTS_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    pack_.ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);

    outputs_ = OutputFuncPtr(new Output());

    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        std::cout << VectorOutputNames_[i] << std::endl;
    }
}

void TSoperation::print()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        outputs_ -> getOutputFuncByName(VectorOutputNames_[i])(VectorOutputFileNames_[i]);
    }
}