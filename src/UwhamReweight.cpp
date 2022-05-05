#include "UwhamReweight.h"

namespace ReweightRegistry
{
    registry<UwhamReweight> registerUwhamReweight("UwhamReweight");
}

UwhamReweight::UwhamReweight(const ReweightInput& input)
:Reweight(input)
{
    output_->registerOutputFunc("averages", [this](std::string name) -> void {this -> printAverages(name);});

    // dynamically check the type of Wham
    Uwham_ = dynamic_cast<Uwham*>(wham_);
    ASSERT((Uwham_ != nullptr), "The wham inputted with type " << wham_ -> type() << \
    " cannot be applied to this reweight strategy.");
}

void UwhamReweight::calculate()
{
    const auto& xi = Uwham_->getxi();
    const auto& lnwji = Uwham_->getlnwji();
    const auto& MapBinToIndex = Uwham_->getMapBinIndexTolnwjiIndex();
    const auto& DataBinIndex = Uwham_->getDataBinIndex();

    // resize the variables
    std::vector<Real> ones(xi.size(),1.0);

    // get 
    dimension_ = xi[0].size();
    averages_.resize(numBias_,std::vector<Real>(dimension_,0.0));
    FE_.resize(numBias_);

    ASSERT((DataBinIndex.size() == xi.size()), "The size of the bindata mismatches the xi data.");

    // Iterate through each of the Biases for reweighting 
    for (int i=0;i<numBias_;i++)
    {
        std::vector<Real> lnpji(xi.size(),0.0);

        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            lnpji[j] = lnwji[j] - value;
        }

        // calculate the normalization constant 
        Real fk = -1.0 * WhamTools::LogSumExpOMP(lnpji_vec, ones);

        #pragma omp parallel
        {
            // normalized lnpji
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                lnpji[j] = lnpji_vec[j] + fk;
            }

            // declare local average
            std::vector<Real> LocalAverage(dimension_,0.0);
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                for (int k=0;k<dimension_;k++)
                {
                    LocalAverage[k] += std::exp(lnpji[j]) * xi[j][k];
                }
            }

            #pragma omp critical
            for (int k=0;k<dimension_;k++)
            {
                averages_[i][k] += LocalAverage[k];
            }
        }
        std::cout << "Done with bias " << i << std::endl;
    }
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