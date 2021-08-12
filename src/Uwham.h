#include "Wham.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>

class Uwham:public Wham
{
    public:
        Uwham(const ParameterPack& pack);
        virtual ~Uwham(){};

        virtual void calculate();
    private:
};