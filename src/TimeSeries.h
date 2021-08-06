#include "tools/InputParser.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>


class TimeSeries
{
    public:
        using Real = CommonTypes::Real;

        TimeSeries(const ParameterPack& pack);
        ~TimeSeries(){};
    
    private:
};