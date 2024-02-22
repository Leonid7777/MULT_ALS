#include <iostream>
#include <complex>
#include <random>
#include <cmath>


void stirling_vector(std::complex<double>* mas, int N, double angle) {

    for (int i = 0; i < N; i++) {
        ((mas)[i]).real(std::cos(angle * i));
        ((mas)[i]).imag(std::sin(angle * i));
    }
}

void
vec_to_conj_vec(std::complex<double>* mas, std::complex<double>* sec, std::complex<double>& sum, int N)
{
    sum = 0;
    for(int i = 0; i < N; i++) {
        sum += mas[i] * (conj(sec[i]));
    }
}

void
find_constant(std::complex<double>* mas, std::complex<double>* res, int N, std::complex<double>& C)
{

    std::complex<double> f_sum = 0, s_sum = 0;
    for(int i = 0; i < N; i++) {
        f_sum += mas[i] * conj(res[i]);
        s_sum += res[i] * conj(res[i]);
    }

    C = f_sum / s_sum;
}

void
vec_app(int h, int N, std::complex<double>* mas, std::complex<double>* res, std::complex<double>& C, int length)
{

    double phi = 2 * M_PI / h;
    int val = 0;
    std::complex<double> max_sum = 0, sum = 0;
    std::complex<double>* sec = new std::complex<double>[N];

    for(int i = 0; i < h; i++) {

        stirling_vector(sec, N, i * phi);

        vec_to_conj_vec(mas, sec, sum, N);

        if(std::abs(sum) > std::abs(max_sum)) {
            max_sum = sum;
            val = i;
        }

    }

    phi *= val;

    stirling_vector(res, length, phi);

    find_constant(mas, res, N, C);

    delete[] sec;
}
