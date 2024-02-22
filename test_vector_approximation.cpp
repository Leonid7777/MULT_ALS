#include <iostream>
#include <complex>
#include <random>
#include <cmath>

extern "C" {
double dznrm2_(const int *, const std::complex<double> *, const int *);
}

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
    int ione = 1;
    int szfull = N;

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

    std::cout << std::endl << C << std::endl;

    for(int i = 0; i < N; i++) {
        sec[i] = mas[i] - res[i];
    }

    std::cout << dznrm2_(&szfull, sec, &ione) / dznrm2_(&szfull, res, &ione) << std::endl;

    delete[] sec;
}

int 
main()
{
    int N = 10;
    int length = 20;
    int h = 20;
    double angle = 97 * M_PI / 200;
    std::complex<double> C = 1;
    double phi = 2 * M_PI / h;
    std::complex<double>* mas = new std::complex<double>[N];
    std::complex<double>* res = new std::complex<double>[length];

    // std::cout << "start_mas ";
    // for(int i = 0; i < N; i++) {
    //     std::cout << mas[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "start_res ";
    // for(int i = 0; i < length; i++) {
    //     std::cout << res[i] << " ";
    // }
    // std::cout << std::endl;


    int val = 0;

    stirling_vector(mas, N, angle);

    vec_app(h, N, mas, res, C, length);

    // std::cout << "end_mas ";
    // for(int i = 0; i < N; i++) {
    //     std::cout << mas[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "end_res ";
    // for(int i = 0; i < length; i++) {
    //     std::cout << res[i] << " ";
    // }
    // std::cout << std::endl;

    std::cout << C << std::endl;




    return 0;
}
