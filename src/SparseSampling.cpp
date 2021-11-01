#include "SparseSampling.h"

namespace WhamRegistry
{
    registry<SparseSampling> registerss("sparsesampling");
}

SparseSampling::SparseSampling(const WhamInput& input)
: Wham(input)
{
    initializeBias();

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        ASSERT((VectorTimeSeries_[i] -> getMean().size() == 1), "The dimension of the data must be 1.");
        means_.push_back(VectorTimeSeries_[i] -> getMean());
        ASSERT((VectorTimeSeries_[i] -> getstd().size() == 1), "The dimension of the data must be 1.");
        std_.push_back(VectorTimeSeries_[i] -> getstd()[0]);
    }

    registerOutput("FE", [this](std::string name) -> void {this -> printFE(name);});
    registerOutput("std", [this](std::string name) -> void {this -> printstd(name);});
    registerOutput("force", [this](std::string name) -> void {this -> printForce(name);});
}

void SparseSampling::calculate()
{
    energy_.resize(Biases_.size(),0.0);
    force_.resize(Biases_.size(),0.0);
    preFactor_.resize(std_.size());
    xstars_.resize(Biases_.size(),0.0);

    for (int i=0;i<Biases_.size();i++)
    {
        energy_[i] = Biases_[i] ->getBeta() * Biases_[i] -> calculate(means_[i]);
        auto f = Biases_[i] -> getBeta() * Biases_[i] -> calculateForce(means_[i])[0];

        // assert size of force is 1
        force_[i] = f;
        ASSERT((Biases_[i] -> getXstar().size() == 1), "The dimension of the bias must be 1.");
        xstars_[i] = Biases_[i] -> getXstar()[0];
    }

    for (int i=0;i<std_.size();i++)
    {
        preFactor_[i] = 0.5 * std::log(2 * PI * std_[i] * std_[i]);
    }

    sparseIntegration_.push_back(0);
    // perform sparse integration 
    for (int i=1;i<VectorTimeSeries_.size();i++)
    {
        Real area_=0.0;
        for (int j=0;j<i;j++)
        {
            area_ += (xstars_[j+1] - xstars_[j]) * (force_[j] + force_[j+1]) * 0.5;
        }
        sparseIntegration_.push_back(area_);
    }

    ASSERT((sparseIntegration_.size() == VectorTimeSeries_.size()), "The sparse sampling integration result size does not match with \
    number of time series.");
    FE_.resize(VectorTimeSeries_.size(),0.0);

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        FE_[i] = preFactor_[i] - energy_[i] + sparseIntegration_[i];
    }

    auto minnum = *std::min_element(FE_.begin(), FE_.end());

    for (int i=0;i<FE_.size();i++)
    {
        FE_[i] -= minnum;
    }
}

void SparseSampling::printstd(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# Number \t std\n";

    for (int i=0;i<std_.size();i++)
    {
        ofs_ << i << "\t" << std_[i] << "\n";
    }

    ofs_.close();
}

void SparseSampling::printForce(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# Number \t force\n";

    for (int i=0;i<force_.size();i++)
    {
        ofs_ << i << "\t" << force_[i] << "\n";
    }

    ofs_.close();
}


void SparseSampling::printFE(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# mean \t prefactor \t energy \t sparseIntegrate \t FE\n";

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        ofs_ << means_[i][0] << "\t" << preFactor_[i] << "\t" << energy_[i] << "\t" << sparseIntegration_[i] << "\t" << FE_[i] << "\n";
    }
    ofs_.close();
}
