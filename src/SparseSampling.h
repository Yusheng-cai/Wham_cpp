#pragma once 

#include "Wham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>

class SparseSampling : public Wham
{
    public:
        SparseSampling(const WhamInput& input);

        // make the assumption that input bias and timeseries has the same order
        virtual void calculate() override;

        void printFE(std::string name);
        void printstd(std::string name);
        void printForce(std::string name);

    private:
        std::vector<std::vector<Real>> means_;
        std::vector<Real> std_;
        std::vector<Real> energy_;
        std::vector<Real> force_;
        std::vector<Real> preFactor_;
        std::vector<Real> sparseIntegration_;
        std::vector<Real> FE_;
        std::vector<Real> xstars_;
        Real PI=3.1415926;
};