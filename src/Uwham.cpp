#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const WhamInput& input)
:Wham(input)
{
    ASSERT((VectorTimeSeries_.size() != 0), "No timeseries data was passed in.");
    N_.resize(VectorTimeSeries_.size());
    std::vector<int> dimensions_(VectorTimeSeries_.size());

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTimeSeries_[i].begin(), VectorTimeSeries_[i].end());
        N_[i] = VectorTimeSeries_[i].getSize(); 
        dimensions_[i] = VectorTimeSeries_[i].getDimension();

        Ntot_ += N_[i];
    }

    for (int i=0;i<dimensions_.size()-1;i++)
    {
        ASSERT((dimensions_[i] == dimensions_[i+1]), "The dimension in the " << i << "th timeseries does not match with the " << i+1 << "th time series");
    }

    // record the dimensions of this Wham calculation
    dimension_ = dimensions_[0];

    auto biases = input.pack_.findParamPacks("bias", ParameterPack::KeyType::Required);
    ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

    for (int i=0;i<biases.size();i++)
    {
        std::string biastype = "simplebias";
        biases[i] -> ReadString("type", ParameterPack::KeyType::Optional, biastype);
        Biasptr b = Biasptr(BiasRegistry::Factory::instance().create(biastype, *biases[i]));

        Biases_.push_back(std::move(b));
    }

    BUki_.resize(Biases_.size(), xi_.size());
    #pragma omp parallel for
    for (int i=0;i<xi_.size();i++)
    {
        for (int j=0;j<Biases_.size();j++)
        {
            Real val = Biases_[j]->calculate(xi_[i]); 
            BUki_(j,i) = Biases_[j]->getBeta()*val;
        }
    }

    // Now read in what type of calculation are you letting Uwham do
    auto whamPack = input.pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    initializeStrat(whamPack);

    bool normRead = whamPack->ReadString("normalizationOutput", ParameterPack::KeyType::Optional, NormalizationFileOutput_);

    if (normRead)
    {
        OpenFile(NormalizationFileofs_, NormalizationFileOutput_);
    }

    auto binPacks = whamPack -> findParamPacks("bins", ParameterPack::KeyType::Required);
    initializeBins(binPacks);

    ASSERT((binPacks.size() == dimension_), "The binning dimension is " << binPacks.size() << " while the dimension of the Wham is " << dimension_);

    bool readpji = whamPack -> ReadString("pjiOutput", ParameterPack::KeyType::Optional, pjiFileOutput_);

    if (readpji)
    {
        OpenFile(pjiFileofs_, pjiFileOutput_);
    }

    bool readlnwji = whamPack -> ReadString("lnwjiOutput", ParameterPack::KeyType::Optional, lnwjiOutput_);
    if (readlnwji)
    {
        OpenFile(lnwjiFileofs_, lnwjiOutput_);
    }

    whamPack -> ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
}

void Uwham::OpenFile(std::ofstream& ofs, std::string& name)
{
    ofs.open(name);

    ASSERT((ofs.is_open()), "The file with name " << name << " is not opened.");
}

void Uwham::initializeBins(const std::vector<const ParameterPack*>& BinPacks)
{
    if (BinPacks.size() != 0)
    {
        for (int i=0;i<BinPacks.size();i++)
        {
            Bins_.push_back(Bin(*BinPacks[i])); 
        }
    }
}


void Uwham::initializeStrat(const ParameterPack* whampack)
{
    std::string strat;
    whampack->ReadString("strategy", ParameterPack::KeyType::Required, strat);

    UwhamStrategyInput input = {BUki_, N_, const_cast<ParameterPack&>(*whampack)};
    strat_ = stratptr(UwhamCalculationStrategyRegistry::Factory::instance().create(strat, input));
}

void Uwham::calculate()
{
    auto start = std::chrono::high_resolution_clock::now();
    strat_ -> calculate();

    auto end = std::chrono::high_resolution_clock::now();

    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Calculation took " << diff.count() << std::endl;

    const auto& lnwji = strat_ -> getlnwji_();

    for (int i=0;i<xi_.size();i++)
    {
        auto& x = xi_[i];
        ASSERT((x.size() == Bins_.size()), "The dimension of the data = " << x.size() << " while the dimension of bins = " << Bins_.size());
        std::vector<int> BinIndex_;
        BinIndex_.resize(Bins_.size());

        bool isInRange = true; 

        for (int j=0;j<Bins_.size();j++)
        {
            // User input is 1 based
            int dim = Bins_[j].getDimension() - 1;

            // If this data point is out of range for one of the bins, it is out of range, so we break
            if (! Bins_[j].isInRange(xi_[i][dim]))
            {
                isInRange = false;
                break;
            }

            int index = Bins_[j].findBin(xi_[i][dim]);
            BinIndex_[j] = index;
        }

        if (isInRange)
        {
            auto it = MapBinIndexToVectorlnwji_.find(BinIndex_);

            if ( it == MapBinIndexToVectorlnwji_.end())
            {
                std::vector<Real> lnwji_vec_;
                lnwji_vec_.push_back(lnwji[i]);
                MapBinIndexToVectorlnwji_.insert(std::make_pair(BinIndex_, lnwji_vec_));
            }
            else
            {
                it -> second.push_back(lnwji[i]);
            }
        }
    }

    for (auto it = MapBinIndexToVectorlnwji_.begin();it != MapBinIndexToVectorlnwji_.end();it++)
    {
        auto& l = it -> second;
        std::vector<Real> ones(l.size(), 1.0);

        Real wji = WhamTools::LogSumExp(l, ones);

        MapBinIndexToWji_.insert(std::make_pair(it->first, wji));
    }

}

void Uwham::printOutput()
{
    if (NormalizationFileofs_.is_open())
    {
        int Nsim = BUki_.getNR();
        NormalizationFileofs_ << std::fixed << std::setprecision(precision_);
        NormalizationFileofs_ << "# normalization constants" << "\n";
        for (int i=0;i<Nsim;i++)
        {
            NormalizationFileofs_ << strat_ -> getFk_()[i] << "\n";
        }

        NormalizationFileofs_.close();
    }

    if (pjiFileofs_.is_open())
    {
        pjiFileofs_ << std::fixed << std::setprecision(precision_);
        pjiFileofs_ << "#";
        for (int i=0;i<Bins_.size();i++)
        {
            pjiFileofs_ << "Position" << Bins_[i].getDimension() << "\t";
        }

        for (int i=0;i<Bins_.size();i++)
        {
            pjiFileofs_ << "Index" << Bins_[i].getDimension() << "\t";
        }
        pjiFileofs_ << "pji\tF";
        pjiFileofs_ << "\n";

        for (auto it = MapBinIndexToWji_.begin();it != MapBinIndexToWji_.end();it++)
        {
            auto& index = it -> first;
            for (int i=0;i<Bins_.size();i++)
            {
                Real pos = Bins_[i].getLocationOfBin(index[i]);
                pjiFileofs_ << pos << " ";
            }

            for (int i=0;i<Bins_.size();i++)
            {
                pjiFileofs_ << index[i] << " ";
            }
            pjiFileofs_ << it -> second << " ";
            pjiFileofs_ << (-1.0)*(it -> second) << "\n"; 
        }
        pjiFileofs_.close();
    }

    if (lnwjiFileofs_.is_open())
    {
        const auto& lnwji = strat_ -> getlnwji_();

        for (int i=0;i<lnwji.size();i++)
        {
            lnwjiFileofs_ << lnwji[i] << std::endl;
        }

        lnwjiFileofs_.close();
    }
}