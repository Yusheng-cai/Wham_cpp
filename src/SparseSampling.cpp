#include "SparseSampling.h"

SparseSampling::SparseSampling(WhamInput& input)
: Wham(input)
{
    initializeBias();

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        ASSERT((VectorTimeSeries_[i] -> getMean().size() == 1), "The dimension of the data must be 1.");
        means_.push_back(VectorTimeSeries_[i] -> getMean());
        std_.push_back(VectorTimeSeries_[i] -> getstd());
    }
}

void SparseSampling::calculate()
{
    energy_.resize(Biases_.size(),0.0);
    force_.resize(Biases_.size(),0.0);

    for (int i=0;i<Biases_.size();i++)
    {
        energy_[i] = Biases_[i] -> calculate(means_[i]);
        auto f = Biases_[i] -> calculateForce(means_[i]);

        ASSERT((f.size() == 1), "The force must be of size 1 right now for sparse sampling.");
        force_[i] = f;
    }

}
