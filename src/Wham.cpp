#include "Wham.h"

#include "Wham.h"

Wham::Wham(const WhamInput& input)
:VectorTimeSeries_(input.VectorTimeSeries_), pack_(input.pack_)
{
    ASSERT((VectorTimeSeries_.size() != 0), "No timeseries data was passed in.");
}

void Wham::registerOutput(std::string name, valueFunction func)
{
    auto it  = MapNameToFunction_.find(name);

    ASSERT((it == MapNameToFunction_.end()), "The output with name " << name << " is already registered.");

    MapNameToFunction_.insert(std::make_pair(name, func));
}

Wham::valueFunction& Wham::printOutputFromName(std::string name)
{
    auto it = MapNameToFunction_.find(name);

    ASSERT(( it != MapNameToFunction_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}

WhamTools::Real WhamTools::LogSumExp(const std::vector<Real>& vector, const std::vector<Real>& N)
{
    ASSERT((vector.size() == N.size()), "The size of the vector is not equal to the size of N.");

    // Find the max of the vector
    auto it = std::max_element(vector.begin(), vector.end());
    Real maxVal = *it;

    Real sum = 0.0;
    for (int i=0;i<vector.size();i++)
    {
        sum += N[i] * std::exp(vector[i]-maxVal);
    }

    sum = std::log(sum) + maxVal;

    return sum;
}

WhamTools::Real WhamTools::LogSumExpOMP(const std::vector<Real>& vector, const std::vector<Real>& N)
{
    ASSERT((vector.size() == N.size()), "The size of the vector is not equal to the size of N.");

    // Find the max of the vector
    auto it = std::max_element(vector.begin(), vector.end());
    Real maxVal = *it;

    Real sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<vector.size();i++)
    {
        sum += N[i] * std::exp(vector[i]-maxVal);
    }

    sum = std::log(sum) + maxVal;

    return sum;
}

WhamTools::Real WhamTools::NormVector(const std::vector<Real>& vector)
{
    Real sum_ = 0.0;

    for (int i=0;i<vector.size();i++)
    {
        sum_ += vector[i]*vector[i];
    }

    return sum_;
}
std::vector<WhamTools::Real> WhamTools::calculatelnWi(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{ 
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    std::vector<Real> lnwji;
    lnwji.resize(Ndata);

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<Ndata;i++)
        {
            std::vector<Real> column;
            column.resize(Nsim);
            for (int j=0;j<Nsim;j++)
            {
                column[j] = fk[j]-1.0*BUki(j,i);
            }

            Real val = WhamTools::LogSumExp(column, N);
            lnwji[i] = -1.0*val;
        }
    }

    return lnwji;
}

std::vector<WhamTools::Real> WhamTools::Gradient(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }

    int Nsim = BUki.getNR();
    int Ndata= BUki.getNC();

    std::vector<Real> ones_(Ndata);
    std::fill(ones_.begin(), ones_.end(), 1.0);

    std::vector<Real> gradient(Nsim);
    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N); 

    std::vector<std::vector<Real>> lnpki(Nsim, std::vector<Real>(Ndata));

    #pragma omp parallel for collapse(2)
    for (int k=0;k<Nsim;k++)
    {
        for (int j=0;j<Ndata;j++)
        {
            lnpki[k][j] = fk[k] - BUki(k,j) + lnwji[j];
        }
    }

    std::vector<Real> lnpk(Nsim);

    for (int k=0;k<Nsim;k++)
    {
        Real val = WhamTools::LogSumExp(lnpki[k], ones_);
        lnpk[k] = val;
    }

    for (int k=0;k<Nsim;k++)
    {
        gradient[k] = -1.0/Ntot*(N[k] - N[k] * std::exp(lnpk[k]));
    }

    return gradient;
}


Matrix<WhamTools::Real> WhamTools::Hessian(const Matrix<Real>& BUki, const std::vector<Real>& fk, const std::vector<Real>& N)
{
    int Ntot = 0;
    for (int i=0;i<N.size();i++)
    {
        Ntot += N[i];
    }
    int Nsim = BUki.getNR();
    int Ndata = BUki.getNC();

    Matrix<Real> Hessian(Nsim, Nsim);

    std::vector<Real> lnwji = calculatelnWi(BUki, fk, N); 
    

    std::vector<std::vector<Real>> pki(Nsim, std::vector<Real>(Ndata));

    #pragma omp parallel for collapse(2)
    for (int k=0;k<Nsim;k++)
    {
        for (int j=0;j<Ndata;j++)
        {
            Real lnpki;
            lnpki = fk[k] - BUki(k,j) + lnwji[j];
            pki[k][j] = std::exp(lnpki);
        }
    }

    for (int i=0;i<Nsim;i++)
    {
        for (int j=0;j<Nsim;j++)
        {
            if (i == j)
            {
                Real sum = 0.0;
                Real sum_sq = 0.0;
                #pragma omp parallel for reduction(+:sum,sum_sq) 
                for (int k=0;k<Ndata;k++)
                {
                    sum += pki[i][k];
                    sum_sq += pki[i][k] * pki[i][k];
                }

                Hessian(i,j) = -1.0/Ntot*(-N[i]*sum + N[i]*N[i]*sum_sq);
            }
            else
            {
                Real sum = 0.0;
                #pragma omp parallel for reduction(+:sum)
                for (int k=0;k<Ndata;k++)
                {
                    sum += pki[i][k] * pki[j][k];
                }

                Hessian(i,j) = -1.0/Ntot*(sum*N[i]*N[j]);
            }
        }
    }

    return Hessian;
}