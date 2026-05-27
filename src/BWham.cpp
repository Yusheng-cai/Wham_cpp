#include "BWham.h"
#include "BwhamBinning.h"

namespace {
std::vector<const Bin*> collectBins(const std::vector<Bwham::BinPtr>& bins)
{
    std::vector<const Bin*> refs;
    refs.reserve(bins.size());

    for (const auto& bin : bins) {
        refs.push_back(bin.get());
    }

    return refs;
}
} // namespace

namespace WhamRegistry {
registry<Bwham> registerBwham("Bwham");
}

Bwham::Bwham(const WhamInput& input) : Wham(input)
{
    initializeBinnedGrid();
    countDataPerBin();
    initializeWil();
    initializeStrategy();

    registerOutput("lnpl", [this](std::string name) -> void { this->printlnpl(name); });
}

void Bwham::printlnpl(std::string name)
{
    std::ofstream ofs_;

    ofs_.open(name);

    ofs_ << std::fixed << std::setprecision(precision_);

    ofs_ << "#Bin \t lnpl \t FE\n";

    for (int i = 0; i < lnpl_.size(); i++) {
        for (int j = 0; j < centerBins_[i].size(); j++) {
            ofs_ << centerBins_[i][j] << "\t";
        }

        ofs_ << lnpl_[i] << "\t" << -lnpl_[i] << "\n";
    }

    ofs_.close();
}

void Bwham::calculate()
{
    strategy_->calculate();

    lnpl_ = strategy_->getlnpl();
}

void Bwham::initializeBinnedGrid()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    auto binPacks = whamPack->findParamPacks("bins", ParameterPack::KeyType::Required);

    bins_.clear();
    centerBins_.clear();
    binIndexToFlat_.clear();

    for (const auto& binPack : binPacks) {
        bins_.push_back(BinPtr(new Bin(*binPack)));
    }

    BwhamBinning::BinGrid grid = BwhamBinning::buildBinGrid(collectBins(bins_), dimension_);
    totalBins_ = grid.totalBins;
    centerBins_ = std::move(grid.centers);
    binIndexToFlat_ = std::move(grid.indexToFlat);
}

void Bwham::initializeWil()
{
    reducedBias_.resize(Biases_.size(), totalBins_);

    for (int i = 0; i < Biases_.size(); i++) {
        for (int j = 0; j < totalBins_; j++) {
            reducedBias_(i, j) = Biases_[i]->getBeta() * Biases_[i]->calculate(centerBins_[j]);
        }
    }
}

void Bwham::initializeStrategy()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    std::string strategyType;
    whamPack->ReadString("strategy", ParameterPack::KeyType::Required, strategyType);

    BwhamStrategyInput input = {reducedBias_, N_, countsPerBin_,
                                const_cast<ParameterPack&>(*whamPack)};
    strategy_ = StrategyPtr(
        BwhamCalculationStrategyRegistry::Factory::instance().create(strategyType, input));
}

void Bwham::countDataPerBin()
{
    countsPerBin_.assign(totalBins_, 0);
    const auto bins = collectBins(bins_);

    for (const auto& sample : xi_) {
        std::vector<int> binIndex;
        if (BwhamBinning::findBinIndexForSample(bins, sample, binIndex)) {
            auto it = binIndexToFlat_.find(binIndex);
            ASSERT((it != binIndexToFlat_.end()), "The BWHAM bin index is not found.");
            countsPerBin_[it->second] += 1;
        }
    }

#ifdef MY_DEBUG
    std::cout << "Printing out countsPerBin" << std::endl;
    for (int i = 0; i < countsPerBin_.size(); i++) {
        std::cout << countsPerBin_[i] << std::endl;
    }
#endif
}
