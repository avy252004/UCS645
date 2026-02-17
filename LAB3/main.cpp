#include <print>
#include <random>
#include <omp.h>
#include "correlate.hpp"

auto generate_random_number() -> double {
        static thread_local std::mt19937 mt{std::random_device{}()};
        static thread_local std::uniform_real_distribution<double> range{-1.0, 1.0};
        return range(mt);
}

auto main(int argc, char** argv) -> int {
        if(argc != 3) {
                std::println("Usage: ./correlate_matrix <number_of_rows> <number_of_columns>");
                return -1;
        }
        int rows{std::atoi(argv[1])};
        int cols{std::atoi(argv[2])};
        int max_threads{omp_get_max_threads()};
        std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
        std::vector<double> flat_array(rows * cols);

        #pragma omp parallel for collapse(2)
        for(int i = 0; i<rows; ++i) {
                for(int j = 0; j<cols; ++j) {
                        double num = generate_random_number();
                        matrix[i][j] = num;
                        flat_array[i*cols + j] = num;
                }
        }

        std::println("Using No Threading");
        auto seq_time{correlate_matrix_sequential(matrix)};
        std::println("Execution time: {:.2f}ms\n",seq_time.count());

        std::println("Using Threading with 2D Heap Allocated Array\n");
        for(int num_threads{2}; num_threads<=max_threads; num_threads += 2){
                std::println("With {} threads", num_threads);
                omp_set_num_threads(num_threads);
                auto unoptimized_par_time{correlate_matrix_parallel_2d_array(matrix)};
                std::println("Execution time: {:.2f}ms", unoptimized_par_time.count());
                std::println("Speed Up: {:.2f}x\n", seq_time.count()/unoptimized_par_time.count());
        }

        std::println("Using Threading with Flat Array\n");
        for(int num_threads{2}; num_threads<=max_threads; num_threads += 2){
                std::println("With {} threads", num_threads);
                omp_set_num_threads(num_threads);
                auto optimized_par_time{correlate_matrix_parallel_2d_array(matrix)};
                std::println("Execution time: {:.2f}ms", optimized_par_time.count());
                std::println("Speed Up: {:.2f}x\n", seq_time.count()/optimized_par_time.count());
        }
        return 0;
}
