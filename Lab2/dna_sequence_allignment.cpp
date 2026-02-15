#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>

constexpr int SEQUQNCE_LENGTH{static_cast<int>(1e4)};
constexpr int MATCH{3};
constexpr int MISMATCH{-3};
constexpr int GAP{-2};
constexpr char POSSIBLE[4] {'A', 'T', 'G', 'C'};
constexpr int tile_sizes[7] {32, 64, 96, 128, 192, 256, 512};

auto get_random_number() -> int {
        static thread_local std::mt19937 mt{std::random_device{}()};
        static thread_local std::uniform_int_distribution<int> range{0, 3};
        return range(mt);
}

auto generate_random_sequence() -> std::string {
        std::string str(SEQUQNCE_LENGTH, 'A');
        #pragma omp parallel for
        for(int i = 0; i<SEQUQNCE_LENGTH; ++i) {
                str[i] = POSSIBLE[get_random_number()];
        }
        return str;
}

auto process_tile(int block_row, int block_col, const std::string& seqA, const std::string& seqB, int tile_size,
                  std::vector<std::vector<int>>& score_matrix, int& global_max_score) -> void {

    int start_i{block_row * tile_size + 1};
    int start_j{block_col * tile_size + 1};

    int end_i{std::min((int)seqA.length() + 1, start_i + tile_size)};
    int end_j{std::min((int)seqB.length() + 1, start_j + tile_size)};

    int local_max{0};

    for (int i = start_i; i < end_i; ++i) {
        for (int j = start_j; j < end_j; ++j) {

            int score_diag = score_matrix[i-1][j-1] + ((seqA[i-1] == seqB[j-1]) ? MATCH : MISMATCH);
            int score_up   = score_matrix[i-1][j] + GAP;
            int score_left = score_matrix[i][j-1] + GAP;

            int current_score = std::max({0, score_diag, score_up, score_left});

            score_matrix[i][j] = current_score;
            local_max = std::max(local_max, current_score);
        }
    }

    if (local_max > global_max_score) {
        #pragma omp critical
        {
            if (local_max > global_max_score) {
                global_max_score = local_max;
            }
        }
    }
}

auto dna_sequqnce_allignment_sequential(const std::string& seqA, const std::string& seqB) -> std::chrono::duration<double, std::milli> {
        int rows{static_cast<int>(seqA.length() + 1)};
        int cols{static_cast<int>(seqB.length() + 1)};
        std::vector<int> prev_row(cols, 0);
        std::vector<int> curr_row(cols, 0);

        int max_score{0};

        auto start{std::chrono::steady_clock::now()};
        for (int i = 1; i < rows; ++i) {
                for (int j = 1; j < cols; ++j) {

                        int score_diag = prev_row[j-1] +
                                ((seqA[i-1] == seqB[j-1]) ? MATCH : MISMATCH);

                        int score_up = prev_row[j] + GAP;

                        int score_left = curr_row[j-1] + GAP;

                        int current_val = std::max({0, score_diag, score_up, score_left});

                        curr_row[j] = current_val;
                        max_score = std::max(max_score, current_val);
                }

                prev_row = curr_row;
                std::fill(curr_row.begin(), curr_row.end(), 0);
        }
        auto end{std::chrono::steady_clock::now()};
        std::chrono::duration<double, std::milli> ms{end-start};
        std::cout << std::format("Score: {}\n", max_score);
        std::cout << std::format("Execution time: {:.2f}ms\n", ms.count());
        return ms;
}

auto dna_sequence_allignment_parallel(const std::string& seqA, const std::string& seqB, int tile_size) -> std::chrono::duration<double, std::milli> {
        int rows{static_cast<int>(seqA.length() + 1)};
        int cols{static_cast<int>(seqB.length() + 1)};

        std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));

        int n_block_rows{(rows - 1 + tile_size - 1) / tile_size};
        int n_block_cols{(cols - 1 + tile_size - 1) / tile_size};

        std::vector<int> deps(n_block_rows * n_block_cols, 0);

        int dummy_root = 0;

        int global_max_score{0};

        auto start{std::chrono::steady_clock::now()};
        #pragma omp parallel
        {
                #pragma omp single
                {
                        for (int i = 0; i < n_block_rows; ++i) {
                                for (int j = 0; j < n_block_cols; ++j) {

                                        int self_idx = i * n_block_cols + j;
                                        int* self = &deps[self_idx];

                                        int* top = (i == 0) ? &dummy_root : &deps[(i-1) * n_block_cols + j];

                                        int* left = (j == 0) ? &dummy_root : &deps[i * n_block_cols + (j-1)];

                                        #pragma omp task \
                                        depend(in: *top) \
                                        depend(in: *left) \
                                        depend(out: *self)
                                        {
                                                process_tile(i, j, seqA, seqB, tile_size, matrix, global_max_score);
                                        }
                                }
                        }
                }
        }
        auto end{std::chrono::steady_clock::now()};
        std::chrono::duration<double, std::milli> ms{end-start};
        std::cout << std::format("Score: {}\n", global_max_score);
        std::cout << std::format("Execution time: {:.2f}ms\n", ms.count());
        return ms;
}

auto main() -> int {
        std::string seqA{generate_random_sequence()};
        std::string seqB{generate_random_sequence()};

        std::cout << "With No thearding\n";
        auto seq_time{dna_sequqnce_allignment_sequential(seqA, seqB)};
        int max_threads{omp_get_max_threads()};
        std::cout << "\nWavefront optimized\n";
        for(int num_threads{2}; num_threads<max_threads; num_threads += 2) {
                omp_set_num_threads(num_threads);
                std::cout << std::format("With {} threads\n", num_threads);
                for(int  i{0}; i<7;++i) {
                        std::cout << std::format("{} x {} tile size\n", tile_sizes[i], tile_sizes[i]);
                        auto par_time{dna_sequence_allignment_parallel(seqA, seqB, tile_sizes[i])};
                        auto speed_up{seq_time.count()/par_time.count()};
                        std::cout << std::format("Speed Up: {:.2f}x\n", speed_up);
                        std::cout << std::format("Throughput: {:.2f}\n", (SEQUQNCE_LENGTH*SEQUQNCE_LENGTH)/par_time.count());
                        std::cout << std::format("Efficiency: {:.2f}\n\n", speed_up/num_threads);
                }
        }
        return 0;
}