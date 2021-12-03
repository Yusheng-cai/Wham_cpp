#include "Uwham.h"

UwhamReweight::UwhamReweight(UwhamReweightInputPack& pack)
:wham_(pack.uwham), pack_(pack.pack)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputs_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, VectorOutputFiles_);
    pack_.ReadNumber("axis", ParameterPack::KeyType::Optional, axis_);

    registerOutputFunc("averages", [this](std::string name) -> void {this -> printAverages(name);});
    registerOutputFunc("FE", [this](std::string name) -> void{this -> printFE(name);});

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

    const auto& xi = wham_.getxi();
    const auto& MapBinToIndex = wham_.getMapBinIndexTolnwjiIndex();

    // get number of bins along the dimension we want to perform conditional probabilities for 
    if (dimension_ > 1)
    {
        numBins_ = wham_.getNumBinsPerDimension(axis_);
    }
    else
    {
        numBins_ = 1;
    }

    // obtain the dimension of the wham calculation
    dimension_ = wham_.getDimension();

    FE_.resize(numBias_, std::vector<Real>(MapBinToIndex.size(),0.0));
    averages_.resize(numBias_, std::vector<Real>(dimension_,0.0));

    ones_.resize(xi.size(), 1.0);
}

void UwhamReweight::calculate()
{
    const auto& xi = wham_.getxi();
    const auto& lnwji = wham_.getlnwji();
    const auto& MapBinToIndex = wham_.getMapBinIndexTolnwjiIndex();

    // find the BUji 
    for (int i=0;i<numBias_;i++)
    {
        std::cout << "Done with bias " << i << std::endl;
        std::vector<Real> lnpji_vec(xi.size(),0.0);

        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            lnpji_vec[j] = lnwji[j] - value;
        }

        Real fk = -1.0 * WhamTools::LogSumExpOMP(lnpji_vec, ones_);

        #pragma omp parallel
        {
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                lnpji_vec[j] = lnpji_vec[j] + fk;
            }

            std::vector<Real> alocal(dimension_,0.0);
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                for (int k=0;k<dimension_;k++)
                {
                    alocal[k] += std::exp(lnpji_vec[j]) * xi[j][k];
                }
            }

            #pragma omp critical
            for (int k=0;k<dimension_;k++)
            {
                averages_[i][k] += alocal[k];
            }
        }

        // calculate free energy
        FE_[i].resize(MapBinToIndex.size());
        int index = 0;
        for(auto it = MapBinToIndex.begin();it != MapBinToIndex.end();it++)
        {
            auto& indices = it ->second;
            std::vector<Real> tempdata_(indices.size(),0.0);
            std::vector<Real> ones(indices.size(),1.0);
            for (int j=0;j<indices.size();j++)
            {
                int ind = indices[j];
                tempdata_[j] = lnpji_vec[ind];
            }

            Real f = WhamTools::LogSumExp(tempdata_,ones);
            FE_[i][index] = -1.0*f;
            index ++;
        }
    }
}

void UwhamReweight::registerOutputFunc(std::string name, outputfunc func)
{
    auto it = MapNameToOutput_.find(name);

    ASSERT((it == MapNameToOutput_.end()), "The output with name " << name << " is already registered.");

    MapNameToOutput_.insert(std::make_pair(name, func));
}

UwhamReweight::outputfunc& UwhamReweight::getOutputByName(std::string name)
{
    auto it = MapNameToOutput_.find(name);

    ASSERT((it != MapNameToOutput_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}

void UwhamReweight::printOutput()
{
    for (int i=0;i<VectorOutputs_.size();i++)
    {
        getOutputByName(VectorOutputs_[i])(VectorOutputFiles_[i]);
    }
}

void UwhamReweight::printFE(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int FE_num = FE_[0].size();

    ofs << "#";

    for (int i=0;i<numBias_;i++)
    {
        ofs << "bias" << i << "\t";
    }
    ofs << "\n";

    for (int i=0;i<FE_num;i++)
    {
        for (int j=0;j<numBias_;j++)
        {
            ofs << FE_[j][i] << "\t";
        }
        ofs << "\n";
    }

    ofs.close();
}


void UwhamReweight::printAverages(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "#data" << "\t";
    for (int i=0;i<dimension_;i++)
    {
        ofs << "dim" << i << "\t";
    }
    ofs << "\n";

    for (int i=0;i<numBias_;i++)
    {
        ofs << i << "\t";
        for (int j=0;j<dimension_;j++)
        {
            ofs << averages_[i][j] << "\t";
        }
        ofs << "\n";
    }

    ofs.close();
}