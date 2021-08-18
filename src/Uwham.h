#pragma once
#include "Wham.h"
#include "tools/InputParser.h"
#include "Array.h"
#include "Bias.h"

#include <vector>
#include <array>
#include <memory>

class Uwham:public Wham
{
    public:
        using Biasptr = std::unique_ptr<Bias>;

        Uwham(const ParameterPack& pack);
        virtual ~Uwham(){};
        Matrix<Real> Hessian(const Matrix<Real>& BUji, const std::vector<Real>& fi, const std::vector<Real>& N);

        virtual void calculate();

    private:
        Matrix<Real> BUji_;
        std::vector<Real> fi_;
        std::vector<std::vector<Real>> xji_;

        std::vector<Biasptr> Biases_;

        int Ntot_;
};