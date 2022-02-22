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
    }

    histogram_.clear();

    histogram_.resize(VectorTS_.size());

    for (int i=0;i<VectorTS_.size();i++)
    {
        auto& Ts = VectorTS_[i];

        // resize to dimension
        histogram_[i].resize(Ts->getDimension());

        // bins is also synonymous with dimension
        for (int j=0;j<Bins_.size();j++)
        {
            auto& b = Bins_[j];
            int dim = b->getDimension() - 1;
            int size = Ts->getSize();

            histogram_[i][j].resize(b->getNumbins(),0.0);

            for (int k=0;k<size;k++)
            {
                if (b->isInRange((*Ts)[k][dim]))
                {
                    int num = b->findBin((*Ts)[k][dim]);
                    histogram_[i][j][num] += 1;
                }
            }
        }
    }
}

void HistogramOP::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int numTs = histogram_.size();
    int dim   = histogram_[0].size();

    ofs << "#";
    for (int i=0;i<numTs;i++)
    {
        for (int j=0;j<dim;j++)
        {
            ofs << "Ts" << i + 1 << "dim" << j + 1 << "\t";
        }
    }

    ofs << "\n";

    for (int i=0;i<dim;i++)
    {
        int numdata = Bins_[i]->getNumbins();
        for (int j=0;j<numdata;j++)
        {
            for(int k=0;k<numTs;k++)
            {
                ofs << histogram_[k][i][j] << "\t";
            }

            ofs << "\n";
        }
    }

    ofs.close();
}