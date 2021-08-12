#include "Eigen/Dense"
#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "TimeSeries.h"
#include "tools/Constants.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <array>

class Wham
{
    public:
        using Real = CommonTypes::Real;

        Wham(const ParameterPack& pack);
        virtual ~Wham(){};

        virtual void calculate() = 0;
    
    protected:
        std::vector<TimeSeries> VectorTimeSeries_;

        // the temperature should be unique with each Wham performed 
        Real temperature_ = 298.15;
        Real beta_;
};

namespace WhamRegistry
{
    using Base = Wham;


}