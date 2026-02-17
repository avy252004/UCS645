#pragma once
#include <chrono>

auto correlate_matrix_sequential(const std::vector<std::vector<double>>&) -> std::chrono::duration<double, std::milli>;
auto correlate_matrix_parallel_2d_array(const std::vector<std::vector<double>>&) -> std::chrono::duration<double, std::milli>;
auto correlate_matrix_parallel_flat_array(const std::vector<double>&) -> std::chrono::duration<double, std::milli>;
