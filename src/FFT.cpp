#include "FFT.h"

void FFT::fft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    output.clear();
    output.resize(data.size());
    std::cout << "Starting fftw" << std::endl;

    fftw_complex *in;
    fftw_complex *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*100);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*100);

    std::cout << "Done getting data" << std::endl;
    plan = fftw_plan_dft_1d(100,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    fftw_destroy_plan(plan);
    std::cout << "Destroyed plan" << std::endl;

    return;
}

void FFT::ifft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    output.clear();
    output.resize(data.size());

    fftw_complex *in,*out;
    fftw_plan plan;
    in = (fftw_complex*) data.data();
    out = (fftw_complex*) output.data();

    plan = fftw_plan_dft_1d(data.size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    fftw_destroy_plan(plan);

    return;
}