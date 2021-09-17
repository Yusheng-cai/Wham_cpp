#include "Uwham.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "Bias.h"

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <map>

class UwhamReweight
{
    public:
        using Real = CommonTypes::Real;
        using Biasptr = std::unique_ptr<Bias>;

        // inputted pack in the whamPack
        UwhamReweight(ParameterPack& pack);

        // inputs are the lnwji weights and the xi points
        void calculate(const std::vector<Real>& lnwji, const std::vector<std::vector<Real>>& xi, const std::map<std::vector<int>, std::vector<int>>& Map);

        void printOutput();

    private:
        // type of the bias, default to simple bias
        std::string type_ = "simplebias";

        Biasptr bias_;

        std::string outputName_;
        std::ofstream ofs_;
};