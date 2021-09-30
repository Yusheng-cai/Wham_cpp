#include "BWham.h"

namespace WhamRegistry
{
    registry<Bwham> registerBwham("Bwham");
}

Bwham::Bwham(const WhamInput& input)
:Wham(input)
{
    // intialize the timeseries 
    initializeTimeSeries();

    // initialize the bins
    initializeBins();

    // Bin the data set 
    bindata();

    // initialize the Bias
    initializeBias();

    // initialize Wil
    initializeWil();

    // initialize the strategy
    initializeStrategy();

    registerOutput("lnpl", [this](std::string name)->void { this -> printlnpl(name);});
}

void Bwham::printlnpl(std::string name)
{
    std::ofstream ofs_;

    ofs_.open(name);

    ofs_ << std::fixed << std::setprecision(precision_);

    ofs_ << "#Bin \t lnpl \t FE\n";

    for (int i=0;i<lnpl_.size();i++)
    {
        for (int j=0;j<centerBins_[i].size();j++)
        {
            ofs_ << centerBins_[i][j] << "\t";
        }

        ofs_ << lnpl_[i] << "\t" << -lnpl_[i] << "\n";
    }

    ofs_.close();
}

void Bwham::calculate()
{
    strat_ -> calculate();

    lnpl_ = strat_ -> getlnpl();
}

void Bwham::initializeBins()
{
    auto WhamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    // bins are required for Binned wham 
    auto BinPacks = WhamPack->findParamPacks("bins", ParameterPack::KeyType::Required);

    TotalBins_ = 1;
    for (int i=0;i<BinPacks.size();i++)
    {
        Bins_.push_back(Binptr(new Bin(*BinPacks[i])));
        TotalBins_ = TotalBins_ * Bins_[i]->getNumbins();
    }

    ASSERT((Bins_.size() == dimension_), "The dimension of bin is " << Bins_.size() << " while dimension for data is " << dimension_);

    // number of Bins is the dimension
    dimension_ = Bins_.size();

    ASSERT((dimension_ <= 2), "Currently not performing any larger dimensions than 2 while user have supplied " << dimension_);

    if (dimension_ == 2)
    {
        int index = 0;
        for(int i=0;i<Bins_[0]->getNumbins();i++)
        {
            for (int j=0;j<Bins_[1]->getNumbins();j++)
            {
                Real outerData = Bins_[0] -> getLocationOfBin(i);
                Real innerData = Bins_[1] -> getLocationOfBin(j);
                int dim1 = Bins_[0] -> getDimension() - 1;
                int dim2 = Bins_[0] -> getDimension() - 1;

                std::vector<Real> data(2);
                std::vector<int> BinIndex_(2);

                data[dim1] = outerData;
                data[dim2] = innerData;
                BinIndex_[dim1] = i;
                BinIndex_[dim2] = j;

                centerBins_.push_back(data);
                auto it = MapBinIndexToIndex_.find(BinIndex_);

                ASSERT((it == MapBinIndexToIndex_.end()), "There is duplicate in the bin index.");
                MapBinIndexToIndex_.insert(std::make_pair(BinIndex_,index));
                index ++;
            }
        }
    }
    else
    {
        for(int i=0;i<Bins_[0] -> getNumbins();i++)
        {
            std::vector<Real> data(1,Bins_[0] -> getLocationOfBin(i));
            std::vector<int> BinIndex_(1, i);

            centerBins_.push_back(data);
            auto it = MapBinIndexToIndex_.find(BinIndex_);
            ASSERT((it == MapBinIndexToIndex_.end()), "There is duplicate bin index.");

            MapBinIndexToIndex_.insert(std::make_pair(BinIndex_,i));
        }
    }

    // for (int i=0;i<centerBins_.size();i++)
    // {
    //     for (int j=0;j<dimension_;j++)
    //     {
    //         std::cout << centerBins_[i][j] << "\t";
    //     }
    //     std::cout << "\n";
    // }
}

void Bwham::initializeWil()
{
    BWil_.resize(Biases_.size(), TotalBins_);

    for(int i=0;i<Biases_.size();i++)
    {
        for(int j=0;j<TotalBins_;j++)
        {
            BWil_(i,j) = Biases_[i]->getBeta() * Biases_[i] -> calculate(centerBins_[j]);
        }
    }
}

void Bwham::initializeStrategy()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    std::string strategyType;
    whamPack->ReadString("strategy", ParameterPack::KeyType::Required, strategyType);

    BwhamStrategyInput input = {BWil_, N_, Ml_, const_cast<ParameterPack&>(*whamPack)};
    strat_ = stratptr(BwhamCalculationStrategyRegistry::Factory::instance().create(strategyType, input));
}

void Bwham::bindata()
{
    // Resize Ml to be size of total number of bins
    Ml_.resize(TotalBins_,0);

    for (int i=0;i<xi_.size();i++)
    {
        std::vector<int> BinsIndex(Bins_.size());
        bool inRange = true;
        for(int j=0;j<Bins_.size();j++)
        {
            int dim = Bins_[j] -> getDimension() - 1;
            if (Bins_[j]->isInRange(xi_[i][j]))
            {
                BinsIndex[j] = Bins_[j] -> findBin(xi_[i][j]);
            }
            else
            {
                inRange = false;
                break;
            }
        }

        if (inRange)
        {
            auto it = MapBinIndexToIndex_.find(BinsIndex);

            ASSERT((it != MapBinIndexToIndex_.end()), "The bins index is not found");
            int index = it -> second;

            Ml_[index] += 1;
        }
    }

    // check if there is any Ml that are zero

    for (int i=0;i<Ml_.size();i++)
    {
        ASSERT((Ml_[i] >= 0), "The number of data at bin " << i << " is 0.");
    }

    #ifdef MY_DEBUG
    std::cout << "Printing out Ml" << std::endl;
    for (int i=0;i<Ml_.size();i++)
    {
        std::cout << Ml_[i] << std::endl;
    }
    #endif
}