#include "UwhamReweight.h"

UwhamReweight::UwhamReweight(ParameterPack& pack)
{
    // read the bias type
    pack.ReadString("type", ParameterPack::KeyType::Optional,type_);

    // This param pack contains everything that is needed by the bias
    bias_ = Biasptr(BiasRegistry::Factory::instance().create(type_, pack));

    bool outputRead = pack.ReadString("output", ParameterPack::KeyType::Optional, outputName_);
    if (outputRead)
    {
        ofs_.open(outputName_);
    }
}

void UwhamReweight::calculate(const std::vector<Real>& lnwji, const std::vector<std::vector<Real>>& xi, const std::map<std::vector<int>, std::vector<int>>& Map)
{
    // first we reweight the lnwji
    std::vector<Real> lnwji_reweight(lnwji.size(), 0.0);
    Real beta = bias_ -> getBeta();

    for (int i=0;i<lnwji.size();i++)
    {
        lnwji_reweight[i] = lnwji[i] - beta * bias_ -> calculate(xi[i]);
    }


    // iterate through the Map 
    for (auto it = Map.begin(); it != Map.end(); it ++)
    {
        // The second indices is the indices of lnwji or x 
    }

}