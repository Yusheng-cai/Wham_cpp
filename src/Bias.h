#include "tools/Assert.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>

// base class for all the biases
class Bias
{
    public:
        Bias(const ParameterPack& pack);
        virtual ~Bias(){};

        virtual void calculate() = 0;
    
    private:
};