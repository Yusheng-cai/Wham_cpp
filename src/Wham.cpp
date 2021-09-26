#include "Wham.h"

#include "Wham.h"

Wham::Wham(const WhamInput& input)
:VectorTimeSeries_(input.VectorTimeSeries_), pack_(input.pack_)
{
    ASSERT((VectorTimeSeries_.size() != 0), "No timeseries data was passed in.");

    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    whamPack -> ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    whamPack -> ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);
    whamPack -> ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    ASSERT((VectorOutputNames_.size() == VectorOutputFileNames_.size()), "The output and the output files size is different.");
}

void Wham::registerOutput(std::string name, valueFunction func)
{
    auto it  = MapNameToFunction_.find(name);

    ASSERT((it == MapNameToFunction_.end()), "The output with name " << name << " is already registered.");

    MapNameToFunction_.insert(std::make_pair(name, func));
}

void Wham::initializeTimeSeries()
{
    dimensions_.resize(VectorTimeSeries_.size());
    N_.resize(VectorTimeSeries_.size());

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTimeSeries_[i]->begin(), VectorTimeSeries_[i]->end());
        N_[i] = VectorTimeSeries_[i]->getSize(); 
        std::cout << "Length of data for " << i << " is " << N_[i] << std::endl;
        dimensions_[i] = VectorTimeSeries_[i]->getDimension();

        Ntot_ += N_[i];
    }

    for (int i=0;i<dimensions_.size()-1;i++)
    {
        ASSERT((dimensions_[i] == dimensions_[i+1]), "The dimension in the " << i << "th timeseries does not match with the " << i+1 << "th time series");
    }

    // record the dimensions of this Wham calculation
    dimension_ = dimensions_[0];
}

void Wham::initializeBias()
{
    auto biases = pack_.findParamPacks("bias", ParameterPack::KeyType::Required);

    ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

    for (int i=0;i<biases.size();i++)
    {
        std::string biastype = "simplebias";
        biases[i] -> ReadString("type", ParameterPack::KeyType::Optional, biastype);
        Biasptr b = Biasptr(BiasRegistry::Factory::instance().create(biastype, *biases[i]));

        Biases_.push_back(std::move(b));
    }
}

Wham::valueFunction& Wham::printOutputFromName(std::string name)
{
    auto it = MapNameToFunction_.find(name);

    ASSERT(( it != MapNameToFunction_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}

void Wham::printOutput()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }
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

std::vector<WhamTools::Real> WhamTools::calculatelnpl(const Matrix<Real>& BWil, const std::vector<Real>& Ml, const std::vector<Real>& N, \
    const std::vector<Real>& fk)
{
    int Nbins = Ml.size();
    int Nsim  = BWil.getNR();

    std::vector<Real> lnpl(Nbins, 0.0);

    std::vector<Real> temp;
    for (int i=0;i<Nbins;i++)
    {
        std::vector<Real> temp(Nsim,0.0);
        for (int j=0;j<Nsim;j++)
        {
            temp[j] = fk[j] - BWil(j,i);
        }

        Real res = WhamTools::LogSumExp(temp, N);
        lnpl[i] = std::log(Ml[i]) - res;
    }

    return lnpl;
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
std::vector<WhamTools::Real> WhamTools::BGradient(const Matrix<Real>& BWil, const std::vector<Real>& Ml, const std::vector<Real>& N, const std::vector<Real>& fk)
{
    int Nbins = Ml.size();
    int Nsim  = BWil.getNR();

    Real Ntot = 0;
    for (int i=0;i<Nsim;i++)
    {
        Ntot += N[i];
    }

    // calculate lnpl
    std::vector<Real> lnpl = WhamTools::calculatelnpl(BWil, Ml, N, fk); 

    // calculate the derivative 
    std::vector<Real> ones(Nbins,1.0);
    Real valu = WhamTools::LogSumExp(lnpl, ones);

    std::vector<Real> derivative(Nsim, 0.0);
    for (int i=0;i<Nsim;i++)
    {
        std::vector<Real> temp(Nbins,0.0);
        std::vector<Real> ones(Nbins,1.0);
        for (int j=0;j<Nbins;j++)
        {
            temp[j] = lnpl[j] - BWil(i,j);
        }

        Real res = WhamTools::LogSumExp(temp, ones);
        derivative[i] = N[i] * (std::exp(fk[i] + res) - 1);
    }

    return derivative;
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