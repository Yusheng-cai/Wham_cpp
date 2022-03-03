#include "Wham.h"

Wham::Wham(const WhamInput& input)
:VectorTimeSeries_(input.VectorTimeSeries_), pack_(input.pack_)
{
    ASSERT((VectorTimeSeries_.size() != 0), "No timeseries data was passed in.");

    whamPack_ =  const_cast<ParameterPack*>(pack_.findParamPack("wham", ParameterPack::KeyType::Required));
    whamPack_ -> ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    whamPack_ -> ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);
    whamPack_ -> ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    whamPack_ -> ReadString("name" , ParameterPack::KeyType::Optional, name_);

    ASSERT((VectorOutputNames_.size() == VectorOutputFileNames_.size()), "The output and the output files size is different.");
    registerOutput("histogram", [this](std::string name) -> void{this -> printTimeSeriesBins(name);});
    registerOutput("Autocorrelation", [this](std::string name) -> void{this -> printAutocorrelation(name);});
    registerOutput("forces", [this](std::string name) -> void {this -> printForce(name);});
    registerOutput("Averages", [this](std::string name) -> void {this -> printAverage(name);});

    // initialize the biases
    initializeBias();

    // let's first initialize the time series 
    initializeTimeSeries();

    // Now initialize the bins 
    initializeBins();

    // bin the time series 
    binTimeSeries();
}

void Wham::printAverage(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Timeseries ";
    for (int i = 0; i<dimension_;i++)
    {
        ofs << "OP" << i+1 << " ";
    }
    ofs << "\n";

    for (int i=0;i<Averages_.size();i++)
    {
        ofs << i+1 << " ";
        for (int j=0;j<Averages_[i].size();j++)
        {
            ofs << Averages_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void Wham::printAutocorrelation(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    std::vector<std::vector<Real>> LagTimes_;

    for (auto& ts : VectorTimeSeries_)
    {
        ts -> calculateAutoCorrelation();

        LagTimes_.push_back(ts -> getLagTime());
    }

    ofs << "# ";
    for (int i=0;i<dimension_;i++)
    {
        ofs << "dim" << i+1 << " ";
    }
    ofs << "\n";

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        ofs << i+1 << " ";
        for (int j=0;j<dimension_;j++)
        {
            ofs << LagTimes_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void Wham::printTimeSeriesBins(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int numTs = histogram_.size();
    int dim   = histogram_[0].size();

    ofs << "#";
    for (int i=0;i<numTs;i++)
    {
        for (int j=0;j<dim;j++)
        {
            ofs << "Ts" << i + 1 << "dim" << j + 1 << "\t";
        }
    }

    ofs << "\n";

    for (int i=0;i<dim;i++)
    {
        int numdata = Bins_[i].getNumbins();
        for (int j=0;j<numdata;j++)
        {
            for(int k=0;k<numTs;k++)
            {
                ofs << histogram_[k][i][j] << "\t";
            }

            ofs << "\n";
        }
    }



    ofs.close();
}

void Wham::registerOutput(std::string name, valueFunction func)
{
    auto it  = MapNameToFunction_.find(name);

    ASSERT((it == MapNameToFunction_.end()), "The output with name " << name << " is already registered.");

    MapNameToFunction_.insert(std::make_pair(name, func));
}

void Wham::initializeBins()
{
    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    auto BinPacks = whamPack -> findParamPacks("bins", ParameterPack::KeyType::Required);

    ASSERT((BinPacks.size() == dimension_), "The binning dimension is " << BinPacks.size() << " while the dimension of the Wham is " << dimension_);

    if (BinPacks.size() != 0)
    {
        for (int i=0;i<BinPacks.size();i++)
        {
            Bins_.push_back(Bin(*BinPacks[i])); 
        }
    }
}

void Wham::binTimeSeries()
{
    histogram_.clear();

    histogram_.resize(VectorTimeSeries_.size());

    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        auto Ts = VectorTimeSeries_[i];

        // resize to dimension
        histogram_[i].resize(Ts->getDimension());

        // bins is also synonymous with dimension
        for (int j=0;j<Bins_.size();j++)
        {
            int dim = Bins_[j].getDimension() - 1;
            int size = Ts->getSize();
            auto& b = Bins_[j];

            histogram_[i][j].resize(b.getNumbins(),0.0);

            for (int k=0;k<size;k++)
            {
                if (b.isInRange((*Ts)[k][dim]))
                {
                    int num = b.findBin((*Ts)[k][dim]);
                    histogram_[i][j][num] += 1;
                }
            }
        }
    }
}

void Wham::isRegistered()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        // check if vector output names is registered
        auto it = MapNameToFunction_.find(VectorOutputNames_[i]);

        ASSERT((it != MapNameToFunction_.end()), "The output with name " << VectorOutputNames_[i] << " is not registered.");
    }
}

void Wham::initializeTimeSeries()
{
    ASSERT((VectorTimeSeries_.size() == 1 || VectorTimeSeries_.size() == Biases_.size()), "You can either provided a time series with \
    all the data combined or time series that are equal to size of biases while there are " << Biases_.size() << " biases but just \
    " << VectorTimeSeries_.size() << " time series.");

    // combine the time series into xi
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        xi_.insert(xi_.end(),VectorTimeSeries_[i]->begin(), VectorTimeSeries_[i]->end());
        std::cout << "Length of data for " << i << " is " << VectorTimeSeries_[i]->getSize() << std::endl;
    }

    // if vector time series is not passed in as 1
    if (VectorTimeSeries_.size() > 1)
    {
        std::cout << "Performing uncombined data input." << std::endl;
        dimensions_.resize(VectorTimeSeries_.size());
        N_.resize(VectorTimeSeries_.size());

        for (int i=0;i<VectorTimeSeries_.size();i++)
        {
            N_[i] = VectorTimeSeries_[i] -> getSize();
            Ntot_ += N_[i];
            dimensions_[i] = VectorTimeSeries_[i] -> getDimension();
        }
        
        for (int i=0;i<dimensions_.size()-1;i++)
        {
            ASSERT((dimensions_[i] == dimensions_[i+1]), "The dimension in the " << i << "th timeseries does not match with the " << i+1 << "th time series");
        }

        // record the dimensions of this Wham calculation
        dimension_ = dimensions_[0];
    }
    else
    {
        std::cout << "Performing combined data input." << std::endl;
        int N;
        bool readN = whamPack_->ReadNumber("N", ParameterPack::KeyType::Optional, N);

        dimension_ = VectorTimeSeries_[0] -> getDimension();
        dimensions_.resize(Biases_.size(), dimension_);
        N_.resize(Biases_.size(), N);

        bool readNvec = whamPack_->ReadVectorNumber("Nvec", ParameterPack::KeyType::Optional, N_);

        ASSERT((N_.size() == Biases_.size()), "The inputted N size is different from bias size.");
        ASSERT((readN || readNvec), "Must provide a N in combined data input.");

        Ntot_ = std::accumulate(N_.begin(), N_.end(), Ntot_);
        ASSERT((Ntot_ == xi_.size()), "The inputted Nvec or N does not sum up to the size of xi, one is " << Ntot_ << " while latter is " << xi_.size());
    }

    for (auto ts : VectorTimeSeries_)
    {
        Averages_.push_back(ts->getMean());
        Std_.push_back(ts->getstd());
    }
}

void Wham::printForce(std::string name)
{
    ASSERT((VectorTimeSeries_.size() == Biases_.size()) && (Averages_.size() == Biases_.size()), "To perform the force output, you cannot use the combined output option.");

    std::vector<std::vector<Real>> Forces;

    for (int i=0;i<Biases_.size();i++)
    {
        Forces.push_back(Biases_[i] -> calculateForce(Averages_[i]));
    }

    std::ofstream ofs;
    ofs.open(name);

    ofs << "# ";
    for (int i=0;i<dimension_;i++)
    {
        ofs << "Average" << i+1 << " ";
    }

    for (int i=0;i<dimension_;i++)
    {
        ofs << "Std" << i+1 << " ";
    }

    for (int i=0;i<dimension_;i++)
    {
        ofs << "dFOP" << i+1 << " ";
    }
    ofs << "\n";

    for (int i=0;i<Forces.size();i++)
    {
        for (auto a : Averages_[i])
        {
            ofs << a << " ";
        }

        for (auto s : Std_[i])
        {
            ofs << s << " ";
        }

        for (auto num : Forces[i])
        {
            ofs << num  << " ";
        }
        
        ofs << "\n";
    }

    ofs.close();
}

void Wham::initializeBias()
{
    auto biases = pack_.findParamPacks("bias", ParameterPack::KeyType::Required);

    // no longer needed 
    // ASSERT((biases.size() == VectorTimeSeries_.size()), "The number of time series does not match the number of biases.");

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
        //gradient[k] = -1.0/Ntot*(N[k] - N[k] * std::exp(lnpk[k]));
        gradient[k] = -(N[k] - N[k] * std::exp(lnpk[k]));
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

                //Hessian(i,j) = -1.0/Ntot*(-N[i]*sum + N[i]*N[i]*sum_sq);
                Hessian(i,j) = -1.0*(-N[i]*sum + N[i]*N[i]*sum_sq);
            }
            else
            {
                Real sum = 0.0;
                #pragma omp parallel for reduction(+:sum)
                for (int k=0;k<Ndata;k++)
                {
                    sum += pki[i][k] * pki[j][k];
                }

                //Hessian(i,j) = -1.0/Ntot*(sum*N[i]*N[j]);
                Hessian(i,j) = -1.0*(sum*N[i]*N[j]);
            }
        }
    }

    return Hessian;
}

WhamTools::Real WhamTools::CalculateBAR(const std::vector<Real>& w_F, const std::vector<Real>& w_B, Real DeltaF)
{
    Real sizeWF = w_F.size();
    Real sizeWB = w_B.size();

    Real M = std::log(sizeWF/sizeWB);

    //log f(W) = - log [1 + exp((M + W - DeltaF))]
    //          = - log ( exp[+maxarg] [exp[-maxarg] + exp[(M + W - DeltaF) - maxarg]] )
    //          = - maxarg - log(exp[-maxarg] + exp[(M + W - DeltaF) - maxarg])
    //where maxarg = max((M + W - DeltaF), 0) 
    std::vector<Real> logf_F(sizeWF, 0.0);
    std::vector<Real> onesF(sizeWF,1.0);
    #pragma omp parallel for
    for (int i=0;i<(int)sizeWF;i++)
    {
        Real val = M + w_F[i] - DeltaF;
        Real maxarg = std::max(val, 0.0);

        logf_F[i] = -maxarg - std::log(std::exp(-maxarg) + std::exp(val - maxarg));
    }
    Real log_numer = LogSumExpOMP(logf_F, onesF);

    std::vector<Real> logf_B(sizeWB, 0.0);
    std::vector<Real> onesB(sizeWB,1.0);
    #pragma omp parallel for
    for (int i=0;i<(int)sizeWB;i++)
    {
        Real val = - M + w_B[i] + DeltaF;
        Real maxarg = std::max(val, 0.0);

        logf_B[i] = -maxarg - std::log(std::exp(-maxarg) + std::exp(val - maxarg));
    }
    Real log_denom = LogSumExpOMP(logf_B, onesB);

    return log_numer - log_denom;
}

WhamTools::Real WhamTools::EXP(const std::vector<Real>& w_F)
{
    Real size = w_F.size();

    std::vector<Real> ones(size,1.0);
    std::vector<Real> negw_F(size,0.0);

    for (int i=0;i<size;i++)
    {
        negw_F[i] = - w_F[i];
    }

    Real val = LogSumExpOMP(negw_F, ones);
    Real denom = std::log(size);

    return - ( val - denom);
}

WhamTools::Real WhamTools::CalculateDeltaFBarIterative(const std::vector<Real>& w_F, const std::vector<Real>& w_B, int max_iterations, Real tol)
{
    Real DeltaFold = 0.0;
    Real DeltaF = 0.0;
    for (int i=0;i<max_iterations;i++)
    {
        DeltaF = - CalculateBAR(w_F, w_B, DeltaF) + DeltaFold;
        DeltaFold = DeltaF;

        Real relativeChange = std::abs(DeltaFold - DeltaF)/DeltaFold;

        if (relativeChange < tol)
        {
            break;
        }
    }

    return DeltaF;
}

WhamTools::Real WhamTools::CalculateDeltaFBarBisection(const std::vector<Real>& w_F, const std::vector<Real>& w_B, int max_iterations)
{   
    // give an initial guess 
    Real upperB = EXP(w_F);
    Real LowerB = -EXP(w_B);

    Real FUpperB= CalculateBAR(w_F, w_B, upperB);
    Real FLowerB = CalculateBAR(w_F, w_B, LowerB);

    while ((FUpperB * FLowerB) > 0)
    {
        Real average = 0.5 * (upperB + LowerB);
    }
    Real multiple= FUpperB * FLowerB;

    std::cout << "FupperB = " << FUpperB << "\n";
    std::cout << "FlowerB = " << FLowerB << "\n";
    std::cout << "Multiple = " << multiple << "\n";


    ASSERT((multiple<0.0), "The initial guesses must be opposite sign");

    return 1;
}

WhamTools::Real WhamTools::Uwham_NLL_equation(const std::vector<Real>& f_k, const Matrix<Real>& BUki, const std::vector<Real>& N)
{
    int Nsim = BUki.getNR();
    int Ndata= BUki.getNC();

    // Get the total N 
    Real Ntot = VectorOP::VectorSum(N);

    // get the fraction of N/Ntot
    std::vector<Real> N_fraction(Nsim,0);
    for (int i=0;i<Nsim;i++)
    {
        N_fraction[i] = N[i] / Ntot;
    }

    ASSERT((f_k.size() == Nsim), "The dimension of fk does not match that of the number of simulation.");

    // Calculates the first part of the equation
    Real firstPart = 0.0;
    for (int i=0;i<Nsim;i++)
    {
        firstPart += N[i] * f_k[i];
    }

    // Calculates the second part of the equation
    Real secondPart = 0.0;
    #pragma omp parallel for reduction(+:secondPart)
    for (int i=0;i<Ndata;i++)
    {
        std::vector<Real> temp(Nsim);
        for (int j=0;j<Nsim;j++)
        {
            temp[j] = f_k[j] - BUki(j,i);
        }

        secondPart += WhamTools::LogSumExp(temp,N_fraction);
    }

    return -firstPart + secondPart;
}