#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"
#include "Bin.h"
#include "UwhamCalculationStrategy.h"
#include "parallel/OpenMP_buffer.h"
#include "tools/GeneralTemplateTools.h"

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
        void initializeStrat(Matrix<Real>& BUki, std::vector<Real>& N, std::vector<stratptr>& strategies);
        void initializeBUki();
        void bindata(std::vector<std::vector<Real>>& x, std::map<std::vector<int>, std::vector<int>>& map);

        // This needs to be called after calculation of BUki --> calculates the initial guess using BAR method
        void MakeInitialGuess(const Matrix<Real>& BUki, const std::vector<Real>& N, std::vector<Real>& fk);

        // Make indices as to which group each data point belongs to
        void MakeGroupPointMap(const std::vector<Real>& N, std::vector<std::vector<int>>& GroupIndex);

        // calculate free energy from the bins and lnwji
        void calculateFreeEnergy(const std::vector<Real>& lnwji, std::map<std::vector<int>, std::vector<int>>& map, std::map<std::vector<int>, Real>& FE);

        virtual void calculate() override;
        virtual void printOutput() override;
        virtual void finishCalculate() override {};
        virtual std::string type() override {return "Uwham";}

        // output statements
        void printNormalization(std::string name);
        void printPji(std::string name);
        void printlnwji(std::string name);
        void printTimeSeriesBins(std::string name);
        void printderivative(std::string name);
        void printReweightFE(std::string name);
        void printKL(std::string name);
        void printFEdim(std::string name);

        // calculation functions
        void calculateBUki(const std::vector<std::vector<Real>>& xi, Matrix<Real>& BUki);
        // get error
        void calculateError();

        // getters 
        const std::vector<Real>& getlnwji() const {return lnwji_;}
        const std::vector<std::vector<Real>>& getxi() const {return xi_;}
        const std::vector<std::vector<int>>& getBinnedData() const {return binneddata_;}
        const std::map<std::vector<int>, std::vector<int>>& getMapBinIndexTolnwjiIndex() {return MapBinIndexTolnwjiIndex_;}
        int getNumBinsPerDimension(int num);

        // reduce FE to various dimensions
        void ReduceFEDimension();

    private:
        // Beta Uki energy matrix
        Matrix<Real> BUki_;

        // The lnwji that falls within each of the bins
        std::map<std::vector<int>, std::vector<int>> MapBinIndexTolnwjiIndex_;
        std::map<std::vector<int>, Real> MapBinIndexToWji_;

        // specify which bin each of the data falls into 
        std::vector<std::vector<int>> binneddata_;

        // precision of the output
        int precision_=3;

        // the initial guess for f_k
        std::vector<Real> fk_;
        std::vector<std::vector<int>> GroupIndex_;

        // whether or not we are doing BAR initialization
        bool BAR_=false;
        bool Error_=false;
        int ErrorIter_=0;

        // Map from name to strategy
        std::map<std::string, UWhamCalculationStrategy*> MapNameToStrat_;
        std::vector<stratptr> strategies_;
        std::vector<std::string> strategyNames_;

        // lnwji 
        std::vector<Real> lnwji_;
        // lnpji
        std::vector<std::vector<Real>> lnpji_;

        // error vector
        std::map<std::vector<int>, std::vector<Real>> ErrorFEMap_;
        std::map<std::vector<int>, Real> ErrorMap_;
        std::map<std::vector<int>, Real> MeanMap_;

        // reduced Free energy in each of the dimensions
        std::vector<std::map<int, Real>> FE_dim_;

        // The constructed FE from lnpji
        std::vector<std::map<std::vector<int>, Real>> reweightFE_;

        // the KL divergence 
        std::vector<Real> KL_divergence_;
};