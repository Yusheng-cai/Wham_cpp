#include "Reweight.h"

Reweight::Reweight(const ReweightInput& input)
: wham_(input.wham_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Required, outputNames_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Required, outputfileNames_);

    // find the bias packs
    auto biasPacks = pack_.findParamPacks("bias",ParameterPack::KeyType::Required);
    numBias_ = biasPacks.size();

    // This param pack contains everything that is needed by the bias
    for (int i=0;i<biasPacks.size();i++)
    {
        std::string type = "simplebias";
        biasPacks[i] -> ReadString("type", ParameterPack::KeyType::Optional,type);
        Vectorbias_.push_back(Biasptr(BiasRegistry::Factory::instance().create(type, *biasPacks[i])));
    }

    output_ = outputptr(new Output());
}

void Reweight::printOutput()
{
    for (int i=0;i<outputNames_.size();i++)
    {
        output_->getOutputFuncByName(outputNames_[i])(outputfileNames_[i]);
    }
}