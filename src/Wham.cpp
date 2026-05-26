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
    whamPack_ -> Readbool("verbose", ParameterPack::KeyType::Optional, verbose_);

    ASSERT((VectorOutputNames_.size() == VectorOutputFileNames_.size()), "The output and the output files size is different.");
    registerOutput("histogram", [this](std::string name) -> void{this -> printTimeSeriesBins(name);});
    registerOutput("Autocorrelation", [this](std::string name) -> void{this -> printAutocorrelation(name);});
    registerOutput("forces", [this](std::string name) -> void {this -> printForce(name);});
    registerOutput("Averages", [this](std::string name) -> void {this -> printAverage(name);});
    registerOutput("dataFE", [this](std::string name) -> void {this -> printdataFE(name);});

    // initialize the biases
    initializeBias();

    // let's first initialize the time series 
    initializeTimeSeries();

    // Now initialize the bins 
    initializeBins();

    // bin the time series 
    binTimeSeries();
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
    dataFE_.resize(VectorTimeSeries_.size());

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

    // we do the FE for each of the data 
    for (int i=0;i<VectorTimeSeries_.size();i++)
    {
        auto Ts = VectorTimeSeries_[i];

        int TsSize = Ts->getSize();

        for (int j=0;j<TsSize;j++)
        {
            // the index for the bin --> same size as dimension or number of bins 
            std::vector<int> Index(Bins_.size());

            // initially set in range to true
            bool InRange = true;
            for (int k=0;k<Bins_.size();k++)
            {
                auto& b = Bins_[k];
                int dim = b.getDimension()-1;

                if (b.isInRange((*Ts)[j][dim]))
                {
                    int num = b.findBin((*Ts)[j][dim]);
                    Index[dim] = num;
                }
                else
                {
                    InRange=false;
                    break;
                }
            }

            // if data is in range, then we add it to the free energy
            if (InRange)
            {
                auto it  = dataFE_[i].find(Index);
                if (it != dataFE_[i].end())
                {
                    it -> second  += 1.0/TsSize; 
                }
                else
                {
                    dataFE_[i].insert(std::make_pair(Index, 1.0/TsSize));
                }
            }
        }

        for (auto it = dataFE_[i].begin(); it != dataFE_[i].end(); it ++)
        {
            it -> second = -std::log(it -> second);
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
        if (verbose_){
            std::cout << "Length of data for " << i << " is " << VectorTimeSeries_[i]->getSize() << std::endl;
        }
    }

    // if vector time series is not passed in as 1
    if (VectorTimeSeries_.size() > 1)
    {
        if (verbose_){
            std::cout << "Performing uncombined data input." << std::endl;
        }
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
        if (verbose_){
            std::cout << "Performing combined data input." << std::endl;
        }
        combined_input_=true;
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
