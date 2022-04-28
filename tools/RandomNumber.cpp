#include "RandomNumber.h"

void Random::seed()
{
    instance().seedImpl();
}

void Random::seedImpl()
{
    if (user_defined_seed) 
    {
        std::cout << "Seed using user defined seed " << seed_ << std::endl;
        generator_ = std::mt19937(seed_);
    }
    else
    {
        generator_ = std::mt19937(rd());
    }
};


Random::Real Random::DrawUniform()
{
    Real value = instance().distribution_(generator_);
    ASSERT((value > 0.0), "Value must be larger than 0.0");

    return value;
}

Random::Real Random::DrawExponential(Real lambda){
    Real uniform = DrawUniform();
    Real expr    = 0.0;
    if (uniform > 0.0){
        Real q = std::log(1.0/uniform);
        expr   = q/lambda;
    }

    return expr;
}

Random::Real Random::DrawUniform_minmax(Real min, Real max){
    Real uniform = DrawUniform();
    return min + (max - min)*uniform;
}

Random Random::_instance;


std::vector<int> RandomTools::RandomPermute(int N, int num)
{
    std::vector<int> Index(N);
    std::iota(Index.begin(), Index.end(), 0);

    std::random_shuffle(Index.begin(), Index.end());

    std::vector<int> Output(num);
    for (int i=0;i<num;i++)
    {
        Output[i] = Index[i];
    }

    return Output;
}
