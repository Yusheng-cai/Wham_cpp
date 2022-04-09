#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "UwhamCalculationStrategy.h"
#include "parallel/OpenMP_buffer.h"

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

        Uwham(const WhamInput& pack);
        virtual ~Uwham(){};
        void initializeStrat();
        void initializeBUki();

        // This needs to be called after calculation of BUki --> calculates the initial guess using BAR method
        void BARInitialGuess();

        // Make indices as to which group each data point belongs to
        void MakeGroupPointMap();

        virtual void calculate() override;
        virtual void printOutput() override;
        virtual void finishCalculate() override;
        virtual std::string type() override {return "Uwham";}

        // output statements
        void printNormalization(std::string name);
        void printPji(std::string name);
        void printlnwji(std::string name);
        void printTimeSeriesBins(std::string name);
        void printderivative(std::string name);
        void printderivativeNormTS(std::string name);
        void printReweightFE(std::string name);
        void printKL(std::string name);

        // getters 
        const std::vector<Real>& getlnwji() const {return lnwji_;}
        const std::map<std::vector<int>, std::vector<Real>>& getMapBinIndexToVectorlnwji_() const {return MapBinIndexToVectorlnwji_;}
        const std::map<std::vector<int>, std::vector<int>>& getMapBinIndexTolnwjiIndex() const {return MapBinIndexTolnwjiIndex_;}
        const std::vector<std::vector<Real>>& getxi() const {return xi_;}
        const std::vector<std::vector<int>>& getBinnedData() const {return binneddata_;}
        int getNumBinsPerDimension(int num);

    private:
        // Beta Uki energy matrix
        Matrix<Real> BUki_;

        // strategy for solving the Uwham
        stratptr strat_;
        std::vector<stratptr> strategies_;
        std::vector<std::string> strategyNames_;

        // The lnwji that falls within each of the bins
        std::map<std::vector<int>, std::vector<Real>> MapBinIndexToVectorlnwji_;
        std::map<std::vector<int>, std::vector<int>> MapBinIndexTolnwjiIndex_;
        std::map<std::vector<int>, Real> MapBinIndexToWji_;
        OpenMP::OpenMP_buffer<std::map<std::vector<int>, std::vector<Real>>> MapBinIndexToVectorlnwjiBuffer_;
        OpenMP::OpenMP_buffer<std::map<std::vector<int>, std::vector<int>>> MapBinIndexTolnwjiIndexBuffer_;

        // specify which bin each of the data falls into 
        std::vector<std::vector<int>> binneddata_;

        // precision of the output
        int precision_=3;

        // the initial guess for f_k
        std::vector<Real> fk_;
        std::vector<std::vector<int>> GroupIndex_;
        std::vector<Real> cumulativeSumN_;

        // whether or not we are doing BAR initialization
        bool BAR_=false;

        // Map from name to strategy
        std::map<std::string, UWhamCalculationStrategy*> MapNameToStrat_;

        // lnwji 
        std::vector<Real> lnwji_;

        // lnpji
        std::vector<std::vector<Real>> lnpji_;

        // The constructed FE from lnpji
        std::vector<std::map<std::vector<int>, Real>> reweightFE_;

        // the KL divergence 
        std::vector<Real> KL_divergence_;
};