#include <iostream>
#include <complex>
#include <random>
#include <cmath>





void 
stirling_vector(std::complex<double>* mas, int N, double angle) {

    for (int i = 0; i < N; i++) {
        ((mas)[i]).real(std::cos(angle * i));
        ((mas)[i]).imag(std::sin(angle * i));
    }
}

void 
stirling_vector_with_positions(std::complex<double>* mas, int N, double angle, int* positions) {

    for (int i = 0; i < N; i++) {
        ((mas)[i]).real(std::cos(angle * positions[i]));
        ((mas)[i]).imag(std::sin(angle * positions[i]));
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
vec_app(int h, int N, std::complex<double>* mas, std::complex<double>* res, std::complex<double>& C, int length, int* positions)
{

    double phi = 2 * M_PI / h;
    int val = 0;
    std::complex<double> max_sum = 0, sum = 0;
    std::complex<double>* sec = new std::complex<double>[N];

    for(int i = 0; i < h; i++) {

        stirling_vector_with_positions(sec, N, i * phi, positions);

        vec_to_conj_vec(mas, sec, sum, N);

        if(std::abs(sum) > std::abs(max_sum)) {
            max_sum = sum;
            val = i;
        }

    }

    phi *= val;

    std::cout << phi << std::endl;

    stirling_vector(res, length, phi);

    find_constant(mas, res, N, C);

    delete[] sec;
}

void
create_small_stirling_vector(std::complex<double>* mas, std::complex<double>* res, int* positions, int N, int length, double angle)
{
    int num;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrib(0, length - 1);
    std::vector<int> val;
    int count = 0;

    while (count < N) {
        num = distrib(gen);

        if (!std::count(val.begin(), val.end(), num)) {
            val.push_back(num);
            count++;
        }

        std::sort(val.begin(), val.end());
    }

    for(int i = 0; i < N; i++) {
        positions[i] = val[i];
    }

    stirling_vector_with_positions(mas, N, angle, positions);

}

int
main()
{
    int h, N, length;
    h = 200;
    length = 10;
    N = length / 2;
    double angle = 159 * M_PI / 200;
    int* positions = new int[N];
    std::complex<double> C;
    std::complex<double>* mas = new std::complex<double>[N];
    std::complex<double>* res = new std::complex<double>[length];

    create_small_stirling_vector(mas, res, positions, N, length, angle);

    std::cout << angle << std::endl;

    vec_app(h, N, mas, res, C, length, positions);

    for(int i = 0; i < length; i++) {
        std::cout << res[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
