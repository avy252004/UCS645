#include <chrono>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <omp.h>
#include <random>
#include <vector>
using namespace std;

double get_random_number(){
        static thread_local std::mt19937 mt{std::random_device{}()};
        static thread_local std::uniform_real_distribution<double> random_range{0.0,1.0};
        return random_range(mt);
}

int main(int argc, char** argv){
        if(argc < 2){
                cout << "Usage: ./daxpy_loop scalar_value\n";
                return EXIT_FAILURE;
        }
        int SIZE{1<<16};
        int totol_num_threads{omp_get_max_threads()};
        double a{atof(argv[1])};
        vector<double> x(SIZE, 0.0);
        vector<double> y(SIZE, 0.0);

        for(int start_threads{2};start_threads <= totol_num_threads; ++start_threads){
                cout << "Using threads: " << start_threads << '\n';
                omp_set_num_threads(start_threads);
                auto start_time{std::chrono::steady_clock::now()};
                #pragma omp parallel for
                for(int i = 0;i<SIZE;++i){
                        x[i] = get_random_number();
                        y[i] = get_random_number();
                        x[i] = a*x[i] + y[i];
                }
                auto end_time{std::chrono::steady_clock::now()};
                auto exec_time{end_time - start_time};
                chrono::duration<double, std::milli> ms{exec_time};
                cout << "Execution time: " << ms.count() << "ms\n";
        }
        return EXIT_SUCCESS;
}
