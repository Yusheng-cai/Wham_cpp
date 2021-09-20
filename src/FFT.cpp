#include "FFT.h"

void FFT::fft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    int datasize = data.size();
    output.clear();
    output.resize(datasize);

    fftw_complex *in;
    fftw_complex *out;
    fftw_plan plan;

    in = (fftw_complex*) data.data();
    out = (fftw_complex*) output.data();


    plan = fftw_plan_dft_1d(datasize,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    fftw_destroy_plan(plan);

    return;
}

void FFT::ifft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    int datasize = data.size();

    output.clear();
    output.resize(datasize);

    fftw_complex *in,*out;
    fftw_plan plan;

    in = (fftw_complex*) data.data();
    out = (fftw_complex*) output.data();
    
    plan = fftw_plan_dft_1d(datasize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    for (int i=0;i<datasize;i++)
    {
        output[i] /= datasize;
    }

    fftw_destroy_plan(plan);

    return;
}