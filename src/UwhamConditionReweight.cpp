#include "UwhamConditionReweight.h"

namespace ReweightRegistry
{
    registry<UwhamConditionReweight> registerConditionReweight("ConditionalReweight");
}

UwhamConditionReweight::UwhamConditionReweight(const ReweightInput& input)
: Reweight(input)
{
    uwhamptr_ = dynamic_cast<Uwham*>(wham_);

    ASSERT((uwhamptr_ != nullptr), "The reweight procedure does not apply to the wham with name " << uwhamptr_->getName());

    pack_.ReadNumber("axis", ParameterPack::KeyType::Required, axis_);
    pack_.ReadNumber("axis_avg", ParameterPack::KeyType::Required, axisavg_);

    axis_--;
    axisavg_--;

    output_->registerOutputFunc("average", [this](std::string name) -> void {this -> printAverage(name);});

    checkOutputValidity();
}

void UwhamConditionReweight::calculate()
{
    const auto& xi = uwhamptr_->getxi();
    const auto& lnwji = uwhamptr_->getlnwji();
    const auto& bindata = uwhamptr_->getBinnedData();
    int numbins = uwhamptr_->getNumBinsPerDimension(axis_);

    // resize axis data 
    axisdata_.resize(numbins);

    // get the index of data per axis
    for (int i=0;i<bindata.size();i++)
    {
        // this means that the data is not in range
        if (bindata[i].size() != 0)
        {
            int index = bindata[i][axis_];
            axisdata_[index].push_back(i);
        }
    }


    peraxisaverage_.resize(numBias_, std::vector<Real>(numbins, 0.0));
    for (int i=0;i<numBias_;i++)
    {
        // first let's update the new weight
        std::vector<Real> lnpji_vec(xi.size(),0.0);

        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            lnpji_vec[j] = lnwji[j] - value;
        }

        std::vector<Real> ones(lnpji_vec.size(),1.0);
        Real fk = - WhamTools::LogSumExpOMP(lnpji_vec,ones);

        #pragma omp parallel for
        for (int j=0;j<xi.size();j++)
        {
            lnpji_vec[j] = lnpji_vec[j] + fk;
        }

        // let's calculate average per axis 
        for (int j=0;j<axisdata_.size();j++)
        {
            if (axisdata_[j].size() > 0 )
            {
                std::vector<Real> peraxisweight_(axisdata_[j].size());
                std::vector<Real> peraxisdata_(axisdata_[j].size());
                for (int k=0;k<axisdata_[j].size();k++)
                {
                    int index = axisdata_[j][k];
                    peraxisweight_[k] = lnpji_vec[index];
                    peraxisdata_[k] = xi[index][axisavg_];
                }

                std::vector<Real> ones(peraxisweight_.size(),1.0);

                Real fk =  - WhamTools::LogSumExp(peraxisweight_, ones);
                for (int k=0;k<axisdata_[j].size();k++)
                {
                    peraxisweight_[k] = peraxisweight_[k] + fk; 
                }

                Real avg=0.0;
                std::vector<Real> avglogsumexp(peraxisdata_.size());
                for (int k=0;k<peraxisdata_.size();k++)
                {
                    avg += std::exp(peraxisweight_[k]) * peraxisdata_[k];
                }

                peraxisaverage_[i][j] = avg;
            }
        }
        std::cout << "Bias " << i << " is done." << std::endl;
    }
}

void UwhamConditionReweight::printAverage(std::string name)
{
    std::ofstream ofs;

    ofs.open(name);

    for (int i=0;i<peraxisaverage_.size();i++)
    {
        for (int j=0;j<peraxisaverage_[i].size();j++)
        {
            ofs << peraxisaverage_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}