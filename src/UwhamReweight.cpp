#include "UwhamReweight.h"

namespace ReweightRegistry
{
    registry<UwhamReweight> registerUwhamReweight("UwhamReweight");
}

UwhamReweight::UwhamReweight(const ReweightInput& input)
:Reweight(input)
{
    // register the outputs 
    output_->registerOutputFunc("ReweightAverages", [this](std::string name) -> void {this -> printReweightAverages(name);});

    // dynamically check the type of Wham
    Uwham_ = dynamic_cast<Uwham*>(wham_);
    ASSERT((Uwham_ != nullptr), "The wham inputted with type " << wham_ -> type() << \
    " cannot be applied to this reweight strategy.");
}

void UwhamReweight::calculate()
{
    // obtain the necessary data 
    const auto& xi = Uwham_->getxi();
    const auto& lnwji = Uwham_->getlnwji();

    // resize the variables
    std::vector<Real> ones(xi.size(),1.0);

    dimension_ = xi[0].size();
    averages_.resize(numBias_,std::vector<Real>(dimension_,0.0));

    // Iterate through each of the Biases for reweighting 
    for (int i=0;i<numBias_;i++)
    {
        std::vector<Real> lnpji(xi.size(),0.0);

        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++){
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            lnpji[j] = lnwji[j] - value;
        }

        // calculate the normalization constant 
        Real fk = -1.0 * WhamTools::LogSumExpOMP(lnpji, ones);

        #pragma omp parallel
        {
            // normalized lnpji
            #pragma omp for
            for (int j=0;j<xi.size();j++)
            {
                lnpji[j] = lnpji[j] + fk;
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

void UwhamReweight::printReweightAverages(std::string name)
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