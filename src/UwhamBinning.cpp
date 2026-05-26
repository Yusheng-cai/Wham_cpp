#include "Uwham.h"

#include <cmath>

void Uwham::bindata(std::vector<std::vector<Real>>& xi, std::map<std::vector<int>, std::vector<int>>& map, std::vector<std::vector<int>>& DataBinIndex)
{
    // clear the bin index for each of the data point
    DataBinIndex.clear();
    DataBinIndex.resize(xi.size());

    // clear the map
    map.clear();

    OpenMP::OpenMP_buffer<std::map<std::vector<int>, std::vector<int>>> mapBuffer;
    mapBuffer.set_master_object(map);

    #pragma omp parallel
    {
        // access the local lnwji index map
        auto& localmap  = mapBuffer.access_buffer_by_id();

        // clear the local maps
        localmap.clear();

        #pragma omp for
        for (int i=0;i<xi.size();i++){
            std::vector<int> BinIndex;
            BinIndex.resize(Bins_.size());

            bool isInRange = true;

            for (int j=0;j<Bins_.size();j++){
                // User input is 1 based
                int dim = Bins_[j].getDimension() - 1;

                // If this data point is out of range for one of the bins, it is out of range, so we break
                if (! Bins_[j].isInRange(xi[i][dim])){
                    // break from the for loop
                    isInRange = false;
                    break;
                }

                int index = Bins_[j].findBin(xi[i][dim]);
                BinIndex[j] = index;
            }

            if (isInRange){
                // add the indices to the binned data vector
                templatetools::InsertIntoVectorMap(BinIndex, i, localmap);
                DataBinIndex[i] = BinIndex;
            }
        }
    }

    // Perform the same operation but for map bin index to ln index
    for (auto m = mapBuffer.beginworker(); m != mapBuffer.endworker(); m ++)
    {
        for (auto it = m->begin(); it != m -> end(); it++)
        {
            auto Indexit = map.find(it -> first);

            // if the original map also has the bin index
            if (Indexit != map.end())
            {
                Indexit -> second.insert(Indexit -> second.end(), it -> second.begin(), it ->second.end());
            }
            else
            {
                map.insert(std::make_pair(it -> first, it -> second));
            }
        }
    }
}

void Uwham::calculateFreeEnergy(const std::vector<Real>& lnwji, const std::map<std::vector<int>, std::vector<int>>& map, std::map<std::vector<int>, Real>& FE)
{
    FE.clear();

    for (auto it = map.begin();it != map.end();it++)
    {
        auto& l = it -> second;
        std::vector<Real> ones(l.size(), 1.0);
        std::vector<Real> lnwji_bin(l.size(),0.0);

        for (int i=0;i<l.size();i++)
        {
            lnwji_bin[i] = lnwji[l[i]];
        }

        Real wji = WhamTools::LogSumExp(lnwji_bin, ones);

        FE.insert(std::make_pair(it->first, wji));
    }

    return;
}

void Uwham::ReduceFEDimension()
{
    FE_dim_.clear();
    FE_dim_.resize(Bins_.size());

    for (int i=0;i<Bins_.size();i++){
        // map each of the bin num in various dimension to lnwji index
        std::map<int, std::vector<int>> MapDimBinNumTolnwjiIndex;
        for (auto it = MapBinIndexTolnwjiIndex_.begin(); it != MapBinIndexTolnwjiIndex_.end();it++)
        {
            int binnum = it -> first[i];
            auto value = it -> second;

            auto itt = MapDimBinNumTolnwjiIndex.find(binnum);

            if (itt != MapDimBinNumTolnwjiIndex.end())
            {
                itt -> second.insert(itt -> second.end(), value.begin(), value.end());
            }
            else
            {
                MapDimBinNumTolnwjiIndex.insert(std::make_pair(binnum, value));
            }
        }
        // Now using those index, find the -log(sum(exp(lnwji)))
        for (auto it = MapDimBinNumTolnwjiIndex.begin(); it != MapDimBinNumTolnwjiIndex.end();it++)
        {
            std::vector<Real> lnwji_dim;

            for (int j=0;j<it -> second.size();j++)
            {
                lnwji_dim.push_back(lnwji_[it ->second[j]]);
            }

            std::vector<Real> ones(lnwji_dim.size(),1.0);
            Real Fe = -1.0 * WhamTools::LogSumExpOMP(lnwji_dim, ones);

            FE_dim_[i].insert(std::make_pair(it -> first, Fe));
        }
    }
}

int Uwham::getNumBinsPerDimension(int num)
{
    ASSERT((num < Bins_.size()), "The dimension provided is larger than the total number of dimensions which is " << Bins_.size());

    return Bins_[num].getNumbins();
}
