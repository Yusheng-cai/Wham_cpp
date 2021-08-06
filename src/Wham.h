#include "Eigen/Dense"
#include "tools/InputParser.h"
#include "tools/Assert.h"

#include <vector>
#include <array>

class Wham
{
    public:
        Wham(const ParameterPack& pack);
        virtual ~Wham(){};

        virtual void calculate() = 0;
    
    protected:
};