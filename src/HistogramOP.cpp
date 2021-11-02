#include "HistogramOP.h"

namespace timeseriesOP
{
    registry<HistogramOP> registerHist("histogram");
}

HistogramOP::HistogramOP(const TSInput& input)
: TSoperation(input)
{
    auto binPacks = pack_.findParamPacks("bins", ParameterPack::KeyType::Required);

    for (int i=0;i<binPacks.size();i++)
    {
        Bins_.push_back(Binptr(new Bin(*binPacks[i])));
    }
    binDimension_ = Bins_.size();
    outputs_->registerOutputFunc("histogram", [this](std::string name) -> void {this -> printHistogram(name);});
}

void HistogramOP::calculate()
{
    for (int i=0;i<VectorTS_.size();i++)
    {
        ASSERT(VectorTS_[i]->getDimension() == binDimension_, \
        "The dimension of the timeseries must match that of the bins while TS = " \
         << VectorTS_[i]->getDimension() << " but bin = " << binDimension_);

        xi_.insert(xi_.end(), VectorTS_[i]->begin(), VectorTS_[i]->end());
    }


    // now start the binning
    for (int i=0;i<xi_.size();i++)
    {
        bool isIn = true;
        std::vector<int> binNum_;
        for (int j=0;j<Bins_.size();j++)
        {
            if (! Bins_[j] -> isInRange(xi_[i][j]))
            {
                isIn = false;
                break;
            }
            else
            {
                binNum_.push_back(Bins_[j] -> findBin(xi_[i][j]));
            }
        }

        if (isIn)
        {
            auto it = histogramMap_.find(binNum_);

            if (it == histogramMap_.end())
            {
                histogramMap_.insert(std::make_pair(binNum_, 1));
            }
            else
            {
                it -> second += 1;
            }
        }
    }
}

void HistogramOP::printHistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "#";
    for (int i=0;i<binDimension_;i++)
    {
        ofs_ << "dim" << i+1 << " ";
    }
    ofs_ << "\n";

    for (auto it = histogramMap_.begin(); it != histogramMap_.end(); it++)
    {
        for (int i=0;i<it -> first.size();i++)
        {
            ofs_ << it -> first[i] << " ";
        }

        ofs_ << it -> second << "\n";
    }

    ofs_.close();
}