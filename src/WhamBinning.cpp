#include "Wham.h"

#include <cmath>

void Wham::initializeBins()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    auto BinPacks = whamPack -> findParamPacks("bins", ParameterPack::KeyType::Required);

    ASSERT((BinPacks.size() == dimension_), "The binning dimension is " << BinPacks.size() << " while the dimension of the Wham is " << dimension_);

    if (BinPacks.size() != 0)
    {
        for (int i=0;i<BinPacks.size();i++)
        {
            Bins_.push_back(Bin(*BinPacks[i]));
        }
    }
}

void Wham::binTimeSeries()
{
    histogram_.clear();

    histogram_.resize(VectorTimeSeries_.size());
    dataFE_.resize(VectorTimeSeries_.size());

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        auto Ts = VectorTimeSeries_[i];

        // resize to dimension
        histogram_[i].resize(Ts->getDimension());

        // bins is also synonymous with dimension
        for (int j=0;j<Bins_.size();j++)
        {
            int dim = Bins_[j].getDimension() - 1;
            int size = Ts->getSize();
            auto& b = Bins_[j];

            histogram_[i][j].resize(b.getNumbins(),0.0);

            for (int k=0;k<size;k++)
            {
                if (b.isInRange((*Ts)[k][dim]))
                {
                    int num = b.findBin((*Ts)[k][dim]);
                    histogram_[i][j][num] += 1;
                }
            }
        }
    }

    // we do the FE for each of the data
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        auto Ts = VectorTimeSeries_[i];

        int TsSize = Ts->getSize();

        for (int j=0;j<TsSize;j++)
        {
            // the index for the bin --> same size as dimension or number of bins
            std::vector<int> Index(Bins_.size());

            // initially set in range to true
            bool InRange = true;
            for (int k=0;k<Bins_.size();k++)
            {
                auto& b = Bins_[k];
                int dim = b.getDimension()-1;

                if (b.isInRange((*Ts)[j][dim]))
                {
                    int num = b.findBin((*Ts)[j][dim]);
                    Index[dim] = num;
                }
                else
                {
                    InRange=false;
                    break;
                }
            }

            // if data is in range, then we add it to the free energy
            if (InRange)
            {
                auto it  = dataFE_[i].find(Index);
                if (it != dataFE_[i].end())
                {
                    it -> second  += 1.0/TsSize;
                }
                else
                {
                    dataFE_[i].insert(std::make_pair(Index, 1.0/TsSize));
                }
            }
        }

        for (auto it = dataFE_[i].begin(); it != dataFE_[i].end(); it ++)
        {
            it -> second = -std::log(it -> second);
        }
    }
}
