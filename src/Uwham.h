#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "UwhamCalculationStrategy.h"
#include "UwhamReweight.h"

#include <vector>
#include <iomanip>
#include <map>
#include <array>
#include <memory>
#include <iostream>
#include <chrono>

class Uwham:public Wham
{
    public:
        using Biasptr = std::unique_ptr<Bias>;
        using stratptr= std::unique_ptr<UWhamCalculationStrategy>;
        using reweightptr = std::unique_ptr<UwhamReweight>;

        Uwham(const WhamInput& pack);
        virtual ~Uwham(){};
        void initializeStrat();
        void initializeBins();
        void initializeBUki();
        void initializePostProcessing();

        virtual void calculate() override;
        virtual void printOutput() override;
        virtual void finishCalculate() override;

        // output statements
        void printNormalization(std::string name);
        void printPji(std::string name);
        void printlnwji(std::string name);
        void printTimeSeriesBins(std::string name);

        // bin the timeseries 
        void binTimeSeries();

        // getters 
        const UWhamCalculationStrategy& getStrategy() const {return *strat_;}
        const std::vector<Real>& getlnwji() const {return getStrategy().getlnwji_();}
        const std::map<std::vector<int>, std::vector<Real>>& getMapBinIndexToVectorlnwji_() const {return MapBinIndexToVectorlnwji_;}
        const std::map<std::vector<int>, std::vector<int>>& getMapBinIndexTolnwjiIndex() const {return MapBinIndexTolnwjiIndex_;}
        const std::vector<std::vector<Real>>& getxi() const {return xi_;}

    private:
        // Beta Uki energy matrix
        Matrix<Real> BUki_;

        // strategy for solving the Uwham
        stratptr strat_;

        std::vector<Bin> Bins_;

        // The lnwji that falls within each of the bins
        std::map<std::vector<int>, std::vector<Real>> MapBinIndexToVectorlnwji_;
        std::map<std::vector<int>, std::vector<int>> MapBinIndexTolnwjiIndex_;
        std::map<std::vector<int>, Real> MapBinIndexToWji_;

        // precision of the output
        int precision_=3;

        // histogram for each dimension of data
        std::vector<std::vector<std::vector<Real>>> histogram_;

        // declare a vector of reweight objects 
        std::vector<reweightptr> reweight_;
};