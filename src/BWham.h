#pragma once
#include "Wham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "Array.h"
#include "BwhamCalculationStrategy.h"

#include <vector>
#include <map>
#include <memory>
#include <iomanip>

class Bwham : public Wham
{
    public:
        using Binptr = std::unique_ptr<Bin>; 
        using stratptr = std::unique_ptr<BWhamCalculationStrategy>;

        Bwham(const WhamInput& input);

        virtual void calculate() override;
        virtual std::string type() override {return "Bwham";}

        void printlnpl(std::string name);
    
    private:
        void initializeBinnedGrid();
        void countDataPerBin();
        void initializeWil();
        void initializeStrategy();

        std::vector<Binptr> Bins_;
        std::vector<std::vector<Real>> centerBins_;

        // total number of bins
        int TotalBins_ = 0;

        Matrix<Real> BWil_;

        // map from index of bins to index in centerBins_
        std::map<std::vector<int>, int> MapBinIndexToIndex_;

        // Number of data per bin 
        std::vector<Real> Ml_;

        // pointer to the strategy
        stratptr strat_;

        // store the outputs of the normalization as well as lnpl
        std::vector<Real> lnpl_;
};
