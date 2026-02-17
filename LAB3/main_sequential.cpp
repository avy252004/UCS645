#include "correlate.hpp"
#include <print>
#include <random>
#include <omp.h>

auto generate_random_number() -> double {
        static thread_local std::mt19937 mt{std::random_device{}()};
        static thread_local std::uniform_real_distribution<double> range{-1.0, 1.0};
        return range(mt);
}

auto main(int argc, char** argv) -> int {
        if(argc != 3) {
                std::println("Usage: ./correlate_matrix_sequential <number_of_rows> <number_of_columns>");
                return -1;
        }
        int rows{std::atoi(argv[1])};
        int cols{std::atoi(argv[2])};
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
        auto seq_time{correlate_matrix_sequential(matrix)};
        std::println("Execution time: {:.2f}ms",seq_time.count());
}
