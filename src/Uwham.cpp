#include "Uwham.h"

namespace WhamRegistry
{
    registry<Uwham> registerUwham("Uwham");
}

Uwham::Uwham(const WhamInput& input)
:Wham(input)
{
    N_.resize(VectorTimeSeries_.size());
    registerOutput("normalization", [this](std::string name)-> void {this->printNormalization(name);});
    registerOutput("pji", [this](std::string name)->void{this->printPji(name);});
    registerOutput("lnwji", [this](std::string name)->void{this->printlnwji(name);});
    registerOutput("histogram", [this](std::string name) -> void{this -> printTimeSeriesBins(name);});

    // Read in timeseries
    initializeTimeSeries();

    // read in the binning information
    initializeBins();

    // bin Timeseries
    binTimeSeries();

    // initialize bias and calculate BUki
    initializeBias();
    initializeBUki();
 
    // Now read in what type of calculation are you letting Uwham do
    initializeStrat();

    // initialize the Post Processing of Wham
    initializePostProcessing();

    auto whamPack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);
    whamPack -> ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    whamPack -> ReadVectorString("outputs", ParameterPack::KeyType::Optional, VectorOutputNames_);
    whamPack -> ReadVectorString("outputFile", ParameterPack::KeyType::Optional, VectorOutputFileNames_);

    ASSERT((VectorOutputNames_.size() == VectorOutputFileNames_.size()), "The output and the output files size is different.");
}

void Uwham::initializePostProcessing()
{
    auto reweightPack = pack_.findParamPacks("Reweight", ParameterPack::KeyType::Optional);

    for (int i=0;i<reweightPack.size();i++)
    {
        UwhamReweightInputPack input = {*this,const_cast<ParameterPack&>(*reweightPack[i])};
        reweight_.push_back(reweightptr(new UwhamReweight(input)));
    }
}

void Uwham::initializeBUki()
{
    BUki_.resize(Biases_.size(), xi_.size());
    #pragma omp parallel for
    for (int i=0;i<xi_.size();i++)
    {
        for (int j=0;j<Biases_.size();j++)
        {
            Real val = Biases_[j]->calculate(xi_[i]); 
            BUki_(j,i) = Biases_[j]->getBeta()*val;
        }
    }
}

void Uwham::initializeTimeSeries()
{
    dimensions_.resize(VectorTimeSeries_.size());
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

void Uwham::initializeBias()
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

void Uwham::printTimeSeriesBins(std::string name)
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

void Uwham::binTimeSeries()
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


void Uwham::initializeBins()
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

void Uwham::initializeStrat()
{
    auto whampack = pack_.findParamPack("wham", ParameterPack::KeyType::Required);

    std::string strat;
    whampack->ReadString("strategy", ParameterPack::KeyType::Required, strat);

    UwhamStrategyInput input = {BUki_, N_, const_cast<ParameterPack&>(*whampack)};
    strat_ = stratptr(UwhamCalculationStrategyRegistry::Factory::instance().create(strat, input));
}

void Uwham::calculate()
{
    auto start = std::chrono::high_resolution_clock::now();
    strat_ -> calculate();

    auto end = std::chrono::high_resolution_clock::now();

    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Calculation took " << diff.count() << std::endl;

    const auto& lnwji = strat_ -> getlnwji_(); 

    for (int i=0;i<xi_.size();i++)
    {
        auto& x = xi_[i];
        ASSERT((x.size() == Bins_.size()), "The dimension of the data = " << x.size() << " while the dimension of bins = " << Bins_.size());
        std::vector<int> BinIndex_;
        BinIndex_.resize(Bins_.size());

        bool isInRange = true; 

        for (int j=0;j<Bins_.size();j++)
        {
            // User input is 1 based
            int dim = Bins_[j].getDimension() - 1;

            // If this data point is out of range for one of the bins, it is out of range, so we break
            if (! Bins_[j].isInRange(xi_[i][dim]))
            {
                isInRange = false;
                break;
            }
        
            int index = Bins_[j].findBin(xi_[i][dim]);
            BinIndex_[j] = index;
        }

        if (isInRange)
        {
            auto it = MapBinIndexToVectorlnwji_.find(BinIndex_);
            auto it2= MapBinIndexTolnwjiIndex_.find(BinIndex_);

            if ( it == MapBinIndexToVectorlnwji_.end())
            {
                std::vector<Real> lnwji_vec_;
                std::vector<int> lnwjiIndex_vec_;

                lnwji_vec_.push_back(lnwji[i]);
                lnwjiIndex_vec_.push_back(i);

                MapBinIndexToVectorlnwji_.insert(std::make_pair(BinIndex_, lnwji_vec_));
                MapBinIndexTolnwjiIndex_.insert(std::make_pair(BinIndex_, lnwjiIndex_vec_));
            }
            else
            {
                it -> second.push_back(lnwji[i]);
                it2 -> second.push_back(i);
            }
        }
    }

    for (auto it = MapBinIndexToVectorlnwji_.begin();it != MapBinIndexToVectorlnwji_.end();it++)
    {
        auto& l = it -> second;
        std::vector<Real> ones(l.size(), 1.0);

        Real wji = WhamTools::LogSumExp(l, ones);

        MapBinIndexToWji_.insert(std::make_pair(it->first, wji));
    }
}

void Uwham::printNormalization(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int Nsim = BUki_.getNR();
    ofs << std::fixed << std::setprecision(precision_);
    ofs << "# normalization constants" << "\n";

    for (int i=0;i<Nsim;i++)
    {
        ofs << strat_ -> getFk_()[i] << "\n";
    }
    ofs.close();
}

void Uwham::printlnwji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    const auto& lnwji = strat_ -> getlnwji_();

    for (int i=0;i<lnwji.size();i++)
    {
        ofs << lnwji[i] << std::endl;
    }

    ofs.close();
}

void Uwham::printPji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);
    ofs << "#";

    for (int i=0;i<Bins_.size();i++)
    {
        ofs << "Position" << Bins_[i].getDimension() << "\t";
    }

    for (int i=0;i<Bins_.size();i++)
    {
        ofs << "Index" << Bins_[i].getDimension() << "\t";
    }

    ofs << "pji\tF";
    ofs << "\n";

    for (auto it = MapBinIndexToWji_.begin();it != MapBinIndexToWji_.end();it++)
    {
        auto& index = it -> first;
        for (int i=0;i<Bins_.size();i++)
        {
            Real pos = Bins_[i].getLocationOfBin(index[i]);
            ofs << pos << " ";
        }

        for (int i=0;i<Bins_.size();i++)
        {
            ofs << index[i] << " ";
        }
        ofs << it -> second << " ";
        ofs << (-1.0)*(it -> second) << "\n"; 
    }
    ofs.close();
}

void Uwham::printOutput()
{
    for (int i=0;i<VectorOutputNames_.size();i++)
    {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }

    for (int i =0;i<reweight_.size();i++)
    {
        reweight_[i] -> printOutput();
    }
}

void Uwham::finishCalculate()
{
    for (int i=0;i<reweight_.size();i++)
    {
        reweight_[i] -> calculate();
    }
}