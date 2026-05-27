#pragma once
#include "fftw3.h"
#include "tools/CommonTypes.h"

#include "parallel/OpenMP.h"
#include <array>
#include <complex>
#include <string>
#include <vector>

namespace FFT {
using Real = CommonTypes::Real;
using ComplexReal = std::complex<Real>;

void fft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output);
void ifft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output);
} // namespace FFT
