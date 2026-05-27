#include "Uwham.h"

#include <cmath>
#include <fstream>
#include <iomanip>

void Uwham::printOutput()
{
    for (int i = 0; i < VectorOutputNames_.size(); i++) {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }
}

void Uwham::printErrroFE(std::string name)
{
    if (Error_) {
        for (int i = 0; i < ErrorFE_.size(); i++) {
            std::string fname = StringTools::AppendIndexToFileName(name, std::to_string(i));

            printPji(fname, ErrorFE_[i]);
        }
    }
}

void Uwham::printlnwji(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i = 0; i < lnwji_.size(); i++) {
        ofs << lnwji_[i] << std::endl;
    }

    ofs.close();
}

void Uwham::printPji(std::string name, const std::map<std::vector<int>, Real>& FE)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);
    ofs << "#";

    for (int i = 0; i < Bins_.size(); i++) {
        ofs << "Position" << Bins_[i].getDimension() << "\t";
    }

    for (int i = 0; i < Bins_.size(); i++) {
        ofs << "Index" << Bins_[i].getDimension() << "\t";
    }

    ofs << "pji\tF";

    if (Error_) {
        ofs << "\tStd";
        ofs << "\tMean";
    }

    ofs << "\n";

    for (auto it = FE.begin(); it != FE.end(); it++) {
        auto& index = it->first;
        for (int i = 0; i < Bins_.size(); i++) {
            Real pos = Bins_[i].getLocationOfBin(index[i]);
            ofs << pos << " ";
        }

        for (int i = 0; i < Bins_.size(); i++) {
            ofs << index[i] << " ";
        }

        ofs << it->second << " ";
        ofs << (-1.0) * (it->second);

        if (Error_) {
            auto it = ErrorMap_.find(index);

            if (it != ErrorMap_.end()) {
                ofs << " " << it->second;
            } else {
                ofs << " " << 0;
            }

            auto itM = MeanMap_.find(index);

            if (itM != MeanMap_.end()) {
                ofs << " " << itM->second;
            } else {
                ofs << " " << 0;
            }
        }

        ofs << "\n";
    }
    ofs.close();
}

void Uwham::printPji(std::string name)
{
    printPji(name, FreeEnergy_);
}

void Uwham::printNormalization(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int Nsim = BUki_.getNR();
    ofs << std::fixed << std::setprecision(precision_);
    ofs << "# normalization constants" << "\n";

    for (int i = 0; i < Nsim; i++) {
        ofs << fk_[i] << "\n";
    }
    ofs.close();
}

void Uwham::printKL(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i = 0; i < KL_divergence_.size(); i++) {
        ofs << KL_divergence_[i] << "\n";
    }
    ofs.close();
}

void Uwham::printReweightFE(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i = 0; i < reweightFE_.size(); i++) {
        for (auto it = reweightFE_[i].begin(); it != reweightFE_[i].end(); it++) {
            ofs << i + 1 << " ";
            for (int j = 0; j < it->first.size(); j++) {
                ofs << it->first[j] << " ";
            }
            ofs << it->second << "\n";
        }
    }

    ofs.close();
}

void Uwham::printFEdim(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    ofs << "# dimension bin FE \n";
    for (int i = 0; i < FE_dim_.size(); i++) {
        for (auto it = FE_dim_[i].begin(); it != FE_dim_[i].end(); it++) {
            ofs << i + 1 << " " << it->first << " " << it->second << "\n";
        }
    }

    ofs.close();
}

void Uwham::printderivative(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    // we have to calculate the derivative
    int dim = Bins_.size();

    std::vector<std::vector<Real>> dFdN(Biases_.size(), std::vector<Real>(dim, 0.0));
    std::vector<std::vector<Real>> avgs(Biases_.size(), std::vector<Real>(dim, 0.0));

    // each bias has a point
    for (int i = 0; i < Biases_.size(); i++) {
        std::vector<Real> avg(dim, 0.0);
        for (int j = 0; j < xi_.size(); j++) {
            Real factor = fk_[i] - BUki_(i, j) + lnwji_[j];
            factor = std::exp(factor);

            for (int k = 0; k < dim; k++) {
                avg[k] += factor * xi_[j][k];
            }
        }

        avgs[i] = avg;
        std::vector<Real> dF = Biases_[i]->calculateForce(avg);
        ASSERT((dF.size() == dim), "The dimensions don't match.");

        dFdN[i] = dF;
    }

    for (int i = 0; i < Biases_.size(); i++) {
        for (int j = 0; j < dim; j++) {
            ofs << avgs[i][j] << " ";
        }

        for (int j = 0; j < dim; j++) {
            ofs << dFdN[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}
