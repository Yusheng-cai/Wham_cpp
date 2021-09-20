#include "Uwham.h"

UwhamReweight::UwhamReweight(UwhamReweightInputPack& pack)
:wham_(pack.uwham), pack_(pack.pack)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputs_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, VectorOutputFiles_);

    registerOutputFunc("lnpji", [this](std::string name) -> void{this->printlnpji(name);});
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

    dimension_ = xi[0].size();

    BUji_.resize(numBias_,std::vector<Real>(xi.size(),0.0));
    lnpji_.resize(numBias_,std::vector<Real>(xi.size(),0.0));
    Vectorfk_.resize(numBias_,0.0);
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
        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            BUji_[i][j] = value;
        }
    }

    // reweight lnwji to lnpji
    for (int i=0;i<numBias_;i++)
    {
        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = lnwji[j] - BUji_[i][j];

            lnpji_[i][j] = value;
        }
    }

    // calculating normalization factor
    for (int i=0;i<numBias_;i++)
    {
        Real fk = -1.0*WhamTools::LogSumExp(lnpji_[i], ones_);
        Vectorfk_[i] = fk;
    }

    // normalize the lnpji
    for (int i=0;i<numBias_;i++)
    {
        for (int j=0;j<xi.size();j++)
        {
            lnpji_[i][j] = lnpji_[i][j] + Vectorfk_[i];
        }
    }

    // using lnpji, we can calculate the average of the data
    for (int i=0;i<numBias_;i++)
    {
        #pragma omp parallel
        {
            std::vector<Real> alocal(dimension_,0.0);
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                for (int k=0;k<dimension_;k++)
                {
                    alocal[k] += std::exp(lnpji_[i][j]) * xi[j][k];
                }
            }

            #pragma omp critical
            for (int k=0;k<dimension_;k++)
            {
                averages_[i][k] += alocal[k];
            }
        }
    }

    // once it is normalized, then we find its free energy by using the binned info from wham 
    for (int i=0;i<numBias_;i++)
    {
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
                tempdata_[j] = lnpji_[i][ind];
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

void UwhamReweight::printlnpji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    const auto& xi = wham_.getxi();
    int datasize = xi.size();

    ofs << "#";

    for (int i=0;i<numBias_;i++)
    {
        ofs << i << " ";
    }

    ofs << "\n";

    for (int i=0;i<datasize;i++)
    {
        for (int j=0;j<numBias_;j++)
        {
            ofs << lnpji_[j][i] << "\t";
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