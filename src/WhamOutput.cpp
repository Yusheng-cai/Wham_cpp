#include "Wham.h"

#include <fstream>

void Wham::printdataFE(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i = 0; i < dataFE_.size(); i++) {
        for (auto it = dataFE_[i].begin(); it != dataFE_[i].end(); it++) {
            ofs << i + 1 << " ";
            for (int j = 0; j < it->first.size(); j++) {
                ofs << it->first[j] << " ";
            }
            ofs << it->second << "\n";
        }
    }

    ofs.close();
}

void Wham::printAverage(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Timeseries ";
    for (int i = 0; i < dimension_; i++) {
        ofs << "OP" << i + 1 << " ";
    }
    ofs << "\n";

    for (int i = 0; i < Averages_.size(); i++) {
        ofs << i + 1 << " ";
        for (int j = 0; j < Averages_[i].size(); j++) {
            ofs << Averages_[i][j] << " ";
        }

        for (int j = 0; j < Std_[i].size(); j++) {
            ofs << Std_[i][j] << " ";
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

    for (auto& ts : VectorTimeSeries_) {
        ts->calculateAutoCorrelation();

        LagTimes_.push_back(ts->getLagTime());
    }

    ofs << "# ";
    for (int i = 0; i < dimension_; i++) {
        ofs << "dim" << i + 1 << " ";
    }
    ofs << "\n";

    for (int i = 0; i < VectorTimeSeries_.size(); i++) {
        ofs << i + 1 << " ";
        for (int j = 0; j < dimension_; j++) {
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
    int dim = histogram_[0].size();

    ofs << "#";
    for (int i = 0; i < numTs; i++) {
        for (int j = 0; j < dim; j++) {
            ofs << "Ts" << i + 1 << "dim" << j + 1 << "\t";
        }
    }

    ofs << "\n";

    for (int i = 0; i < dim; i++) {
        int numdata = Bins_[i].getNumbins();
        for (int j = 0; j < numdata; j++) {
            for (int k = 0; k < numTs; k++) {
                ofs << histogram_[k][i][j] << "\t";
            }

            ofs << "\n";
        }
    }

    ofs.close();
}

void Wham::printForce(std::string name)
{
    ASSERT((VectorTimeSeries_.size() == Biases_.size() && Averages_.size() == Biases_.size()),
           "To perform the force output, you cannot use the combined output option.");

    std::vector<std::vector<Real>> Forces;

    for (int i = 0; i < Biases_.size(); i++) {
        Forces.push_back(Biases_[i]->calculateForce(Averages_[i]));
    }

    std::ofstream ofs;
    ofs.open(name);

    ofs << "# ";
    for (int i = 0; i < dimension_; i++) {
        ofs << "Average" << i + 1 << " ";
    }

    for (int i = 0; i < dimension_; i++) {
        ofs << "Std" << i + 1 << " ";
    }

    for (int i = 0; i < dimension_; i++) {
        ofs << "dFOP" << i + 1 << " ";
    }
    ofs << "\n";

    for (int i = 0; i < Forces.size(); i++) {
        for (auto a : Averages_[i]) {
            ofs << a << " ";
        }

        for (auto s : Std_[i]) {
            ofs << s << " ";
        }

        for (auto num : Forces[i]) {
            ofs << num << " ";
        }

        ofs << "\n";
    }

    ofs.close();
}

Wham::valueFunction& Wham::printOutputFromName(std::string name)
{
    auto it = MapNameToFunction_.find(name);

    ASSERT((it != MapNameToFunction_.end()),
           "The output with name " << name << " is not registered.");

    return it->second;
}

void Wham::printOutput()
{
    for (int i = 0; i < VectorOutputNames_.size(); i++) {
        std::string name = VectorOutputNames_[i];

        printOutputFromName(name)(VectorOutputFileNames_[i]);
    }
}
