#include "UwhamConditionalReweight.h"

namespace ReweightRegistry
{
    registry<UwhamConditionalReweight> registerConditionalReweight("ConditionalReweight");
}

UwhamConditionalReweight::UwhamConditionalReweight(const ReweightInput& input)
: Reweight(input)
{
    // check if the Wham passed in is Uwham
    uwhamptr_ = dynamic_cast<Uwham*>(wham_);
    ASSERT((uwhamptr_ != nullptr), "The reweight procedure does not apply to the wham with name " << uwhamptr_->getName());

    // axis is the one where we are calculating averages give so <O|axis>
    pack_.ReadNumber("axis", ParameterPack::KeyType::Required, axis_);
    axis_--;

    // register average
    output_->registerOutputFunc("ConditionalAverage", [this](std::string name) -> void {this -> printConditionalAverage(name);});

    checkOutputValidity();
}

void UwhamConditionalReweight::calculate()
{
    // obtain the necessary data 
    const auto& xi = uwhamptr_->getxi();
    const auto& lnwji = uwhamptr_->getlnwji();
    const auto& DataBinIndex = uwhamptr_->getDataBinIndex();
    int numbins = uwhamptr_->getNumBinsPerDimension(axis_);
    dimension_ = xi[0].size();

    // resize axis data 
    ConditionalIndex_.resize(numbins);

    // get the index of data give axis 
    for (int i=0;i<DataBinIndex.size();i++)
    {
        // this means that the data is not in range
        if (DataBinIndex[i].size() != 0)
        {
            int index = DataBinIndex[i][axis_];
            ConditionalIndex_[index].push_back(i);
        }
    }

    // with size (Numbias of bias, numer of bins, data dimension)
    ConditionalAverage_.resize(numBias_, std::vector<std::vector<Real>>(numbins, std::vector<Real>(dimension_,0.0)));

    // iterate over number of bias 
    for (int i=0;i<numBias_;i++)
    {
        // first let's update the new weight based on the bias 
        std::vector<Real> lnpji(xi.size(),0.0);

        // obtain the biased weights 
        #pragma omp parallel for 
        for (int j=0;j<xi.size();j++)
        {
            Real value = Vectorbias_[i]->getBeta()*Vectorbias_[i]->calculate(xi[j]);
            lnpji[j] = lnwji[j] - value;
        }

        // normalized fk 
        std::vector<Real> ones(lnpji.size(),1.0);
        Real fk = - WhamTools::LogSumExpOMP(lnpji,ones);

        // normalize lnpji
        #pragma omp parallel for
        for (int j=0;j<xi.size();j++)
        {
            lnpji[j] = lnpji[j] + fk;
        }

        // let's calculate average give axis 
        for (int j=0;j<ConditionalIndex_.size();j++)
        {
            if (ConditionalIndex_[j].size() > 0 )
            {
                // fill the weight and data give the dimension
                std::vector<Real> ConditionalWeight(ConditionalIndex_[j].size());
                std::vector<std::vector<Real>> ConditionalData(ConditionalIndex_[j].size());
                for (int k=0;k<ConditionalIndex_[j].size();k++)
                {
                    int index = ConditionalIndex_[j][k];

                    // the conditional weight as well as condition data 
                    ConditionalWeight[k] = lnpji[index];
                    ConditionalData[k] = xi[index];
                }

                // Ones 
                std::vector<Real> ones(ConditionalWeight.size(),1.0);

                // normalize per axis weights 
                Real fk =  - WhamTools::LogSumExp(ConditionalWeight, ones);
                for (int k=0;k<ConditionalWeight.size();k++)
                {
                    ConditionalWeight[k] = ConditionalWeight[k] + fk; 
                }

                // find the average for this particular bin
                std::vector<Real> avg(dimension_,0.0);
                for (int k=0;k<ConditionalData.size();k++)
                {
                    for (int m=0;m<dimension_;m++)
                    {
                        avg[m] += std::exp(ConditionalWeight[k]) * ConditionalData[k][m];
                    }
                }

                ConditionalAverage_[i][j] = avg;
            }
        }
        std::cout << "Bias " << i << " is done." << std::endl;
    }
}

void UwhamConditionalReweight::printConditionalAverage(std::string name)
{
    std::ofstream ofs;

    ofs.open(name);

    int AxisNumber = 1;
    pack_.ReadNumber("AxisNumber", ParameterPack::KeyType::Required, AxisNumber);
    AxisNumber--;

    // How should we print out the conditional averages ?
    // TODO : Figure out how to print out higher dimensional conditional averages 
    for (int i=0;i<ConditionalAverage_.size();i++)
    {
        for (int j=0;j<ConditionalAverage_[i].size();j++)
        {
            ofs << ConditionalAverage_[i][j][AxisNumber] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}