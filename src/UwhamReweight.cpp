#include "UwhamReweight.h"

namespace ReweightRegistry
{
    registry<UwhamReweight> registerUwhamReweight("UwhamReweight");
}

UwhamReweight::UwhamReweight(const ReweightInput& input)
:Reweight(input)
{
    output_->registerOutputFunc("averages", [this](std::string name) -> void {this -> printAverages(name);});
    output_->registerOutputFunc("FE", [this](std::string name) -> void{this -> printFE(name);});

    Uwham_ = dynamic_cast<Uwham*>(wham_);

    ASSERT((Uwham_ != nullptr), "The wham inputted with type " << wham_ -> type() << \
    " cannot be applied to this reweight strategy.");
}

void UwhamReweight::calculate()
{
    const auto& xi = Uwham_->getxi();
    const auto& lnwji = Uwham_->getlnwji();
    const auto& MapBinToIndex = Uwham_->getMapBinIndexTolnwjiIndex();
    const auto& bindata = Uwham_->getBinnedData();

    // resize the variables
    ones_.resize(xi.size(),1.0);
    dimension_ = xi[0].size();
    averages_.resize(numBias_,std::vector<Real>(dimension_,0.0));
    FE_.resize(numBias_);

    ASSERT((bindata.size() == xi.size()), "The size of the bindata mismatches the xi data.");

    // find the BUji 
    for (int i=0;i<numBias_;i++)
    {
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
        std::cout << "Done with bias " << i << std::endl;
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