#include "Uwham.h"

#include <cmath>

void Uwham::calculate()
{
    // Start calculation --> first calculate using the entire data set
    for (int i=0;i<strategies_.size();i++){
        UwhamStrategyResult result = strategies_[i] -> calculate(fk_);
        fk_ = result.fk;
        lnwji_ = result.lnwji;
    }

    // calculate the free energy
    calculateFreeEnergy(lnwji_, MapBinIndexTolnwjiIndex_, FreeEnergy_);

    // We should also reweight the data to each of the simulations
    // resize lnpji to the size of 'Number of biases'
    lnpji_.resize(BUki_.getNR(), std::vector<Real>(lnwji_.size(),0.0));
    reweightFE_.resize(BUki_.getNR());

    // Reweight to various biases
    for (int i=0;i<BUki_.getNR();i++){
        #pragma omp parallel for
        for (int j=0;j<lnwji_.size();j++){
            lnpji_[i][j] = fk_[i] - BUki_(i,j) + lnwji_[j];
        }
    }

    // only do this procedure if we are not doing combined input
    if (! combined_input_){
        for (int i=0;i<BUki_.getNR();i++){
            for (auto it = MapBinIndexTolnwjiIndex_.begin(); it != MapBinIndexTolnwjiIndex_.end(); it ++){
                auto& l  = it -> second;
                std::vector<Real> lnpi;

                for (int j=0;j< l.size(); j++)
                {
                    lnpi.push_back(lnpji_[i][l[j]]);
                }

                std::vector<Real> ones(l.size(), 1.0);
                Real pi = -WhamTools::LogSumExp(lnpi, ones);

                reweightFE_[i].insert(std::make_pair(it -> first, pi));
            }
        }

        // Now let's calculate KL divergence
        KL_divergence_.resize(BUki_.getNR(),0.0);
        for (int i=0;i<BUki_.getNR();i++)
        {
            for (auto it = dataFE_[i].begin(); it != dataFE_[i].end(); it ++)
            {
                auto Index = it -> first;

                Real ref_val = it -> second;
                Real prob = std::exp(-ref_val);
                Real val = reweightFE_[i].find(Index) -> second;

                KL_divergence_[i] += prob * (-ref_val + val);
            }
        }
    }

    // get the FE in each of the dimensions
    ReduceFEDimension();

    // if we want to calculate error
    if (Error_){
        calculateError();
    }
}

void Uwham::calculateError()
{
    for (int i=0;i<ErrorIter_;i++)
    {
        std::vector<std::vector<Real>> X;
        std::vector<Real> N;

        for (int i=0;i<VectorTimeSeries_.size();i++)
        {
            auto& ts = VectorTimeSeries_[i];
            std::vector<std::vector<Real>> sample = ts->getIndependentsample();
            X.insert(X.end(), sample.begin(), sample.end());
            N.push_back(sample.size());
        }

        // calculate BUki
        Matrix<Real> BUki;
        calculateBUki(X, BUki);

        // bin the data
        std::map<std::vector<int>,std::vector<int>> map;
        std::vector<std::vector<int>> dbinIndex;
        bindata(X, map, dbinIndex);

        // make initial guess for fk = -log(Qk)
        std::vector<Real> fk_guess;
        std::vector<Real> lnwji;
        MakeInitialGuess(BUki, N, fk_guess);

        // make new strategies
        std::vector<stratptr> strategies;
        std::map<std::vector<int>, Real> FE;
        initializeStrat(BUki, N, strategies);

        for (int i=0;i<strategies.size();i++)
        {
            UwhamStrategyResult result = strategies[i] -> calculate(fk_guess);
            fk_guess = result.fk;
            lnwji = result.lnwji;
        }

        // calculate the Free Energy
        calculateFreeEnergy(lnwji, map, FE);

        // append to the free energy
        ErrorFE_.push_back(FE);

        for (auto it = FE.begin(); it != FE.end(); it ++)
        {
            std::vector<int> copyK = it -> first;
            templatetools::InsertIntoVectorMap(copyK, it -> second, ErrorFEMap_);
        }
    }

    for (auto it = ErrorFEMap_.begin(); it != ErrorFEMap_.end(); it ++)
    {
        Real mean=0.0;
        Real var =0.0;
        for (int i=0;i<it->second.size();i++)
        {
            mean += it -> second[i];
        }

        mean = mean / ErrorIter_;
        for (int i=0;i<it->second.size();i++)
        {
            var += std::pow((it->second[i] - mean),2.0);
        }

        var = var + (ErrorIter_ - it->second.size()) * (mean*mean);
        var = var / ErrorIter_;

        Real std = std::sqrt(var);

        ErrorMap_.insert(std::make_pair(it -> first, std));
        MeanMap_.insert(std::make_pair(it->first, mean));
    }
}
