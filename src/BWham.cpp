#include "BWham.h"
#include "BwhamBinning.h"

namespace
{
std::vector<const Bin*> collectBins(const std::vector<Bwham::Binptr>& bins)
{
    std::vector<const Bin*> refs;
    refs.reserve(bins.size());

    for (const auto& bin : bins)
    {
        refs.push_back(bin.get());
    }

    return refs;
}
}

namespace WhamRegistry
{
    registry<Bwham> registerBwham("Bwham");
}

Bwham::Bwham(const WhamInput& input)
:Wham(input)
{
    initializeBinnedGrid();
    countDataPerBin();
    initializeWil();
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

void Bwham::initializeBinnedGrid()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    auto binPacks = whamPack->findParamPacks("bins", ParameterPack::KeyType::Required);

    Bins_.clear();
    centerBins_.clear();
    MapBinIndexToIndex_.clear();

    for (const auto& binPack : binPacks)
    {
        Bins_.push_back(Binptr(new Bin(*binPack)));
    }

    BwhamBinning::BinGrid grid = BwhamBinning::buildBinGrid(collectBins(Bins_), dimension_);
    TotalBins_ = grid.totalBins;
    centerBins_ = std::move(grid.centers);
    MapBinIndexToIndex_ = std::move(grid.indexToFlat);
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

void Bwham::countDataPerBin()
{
    Ml_.assign(TotalBins_, 0);
    const auto bins = collectBins(Bins_);

    for (const auto& sample : xi_)
    {
        std::vector<int> binIndex;
        if (BwhamBinning::findBinIndexForSample(bins, sample, binIndex))
        {
            auto it = MapBinIndexToIndex_.find(binIndex);
            ASSERT((it != MapBinIndexToIndex_.end()), "The BWHAM bin index is not found.");
            Ml_[it->second] += 1;
        }
    }

    #ifdef MY_DEBUG
    std::cout << "Printing out Ml" << std::endl;
    for (int i=0;i<Ml_.size();i++)
    {
        std::cout << Ml_[i] << std::endl;
    }
    #endif
}
