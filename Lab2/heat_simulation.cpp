#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <string>
#include <format>

constexpr int N{2000};
constexpr int STEPS{500};
constexpr double ALPHA{0.01};
constexpr double DX{1.0};
constexpr double DT{0.1};
constexpr int CHUNK_SIZE{64};

struct SimResult {
        double duration_ms;
        double total_energy;
};

struct ScheduleConfig {
        omp_sched_t type;
        int chunk_size;
        std::string name;
};

auto run_simulation(int num_threads, omp_sched_t sched_type, int chunk_size) -> SimResult {
        std::vector<double> T(N * N, 0.0);
        std::vector<double> T_next(N * N, 0.0);

        int center{N / 2};
        int radius{N / 10};

        #pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                        if ((i - center)*(i - center) + (j - center)*(j - center) < radius*radius) {
                                T[i * N + j] = 100.0;
                        } else {
                                T[i * N + j] = 0.0;
                        }
                }
        }

        double cx {ALPHA * DT / (DX * DX)};
        double cy {ALPHA * DT / (DX * DX)};

        omp_set_num_threads(num_threads);
        omp_set_schedule(sched_type, chunk_size);

        auto start_time{std::chrono::high_resolution_clock::now()};

        for (int step{0}; step < STEPS; ++step) {

        #pragma omp parallel for schedule(runtime)
                for (int i = 1; i < N - 1; ++i) {
                        for (int j = 1; j < N - 1; ++j) {
                                int idx = i * N + j;
                                int up = (i - 1) * N + j;
                                int down = (i + 1) * N + j;
                                int left = i * N + (j - 1);
                                int right = i * N + (j + 1);

                                T_next[idx] = T[idx] +
                                        cx * (T[up] - 2 * T[idx] + T[down]) +
                                        cy * (T[left] - 2 * T[idx] + T[right]);
                        }
                }
                std::swap(T, T_next);
        }

        auto end_time{std::chrono::high_resolution_clock::now()};
        std::chrono::duration<double, std::milli> duration{end_time - start_time};

        double total_energy{0.0};
        #pragma omp parallel for reduction(+:total_energy) num_threads(num_threads)
        for (int i = 0; i < N * N; ++i) {
                total_energy += T[i];
        }

        return { duration.count(), total_energy };
}

auto main() -> int {
        std::cout << "With No Threading\n";
        SimResult seq_res = run_simulation(1, omp_sched_static, 0);

        std::cout << std::format("Execution time: {:.2f}ms\n", seq_res.duration_ms);
        std::cout << std::format("Total Energy: {:.4e}\n", seq_res.total_energy);

        std::vector<ScheduleConfig> strategies {
                {omp_sched_static, 0, "STATIC"},
                {omp_sched_dynamic, CHUNK_SIZE, std::format("DYNAMIC (Chunk {})", CHUNK_SIZE)},
                {omp_sched_guided, CHUNK_SIZE, std::format("GUIDED (Chunk {})", CHUNK_SIZE)}
        };

        int max_threads{omp_get_max_threads()};
        if (max_threads < 2) max_threads = 12;

        for (int t{2}; t <= max_threads; t += 2) {
                std::cout << std::format("\nWith {} threads\n", t);
                for (const auto& strat : strategies) {
                        SimResult par_res{run_simulation(t, strat.type, strat.chunk_size)};
                        double speed_up{seq_res.duration_ms / par_res.duration_ms};

                        double ops = (double)N * N * STEPS;
                        double throughput = ops / (par_res.duration_ms / 1000.0) / 1e6;

                        double efficiency{speed_up / t};

                        std::cout << std::format("[ {:<20} ]\n", strat.name);
                        std::cout << std::format("Execution time: {:.2f}ms\n", par_res.duration_ms);
                        std::cout << std::format("Speed Up:{:.2f}x\n", speed_up);
                        std::cout << std::format("Throughput: {:.2f}\n", throughput);
                        std::cout << std::format("Efficiency: {:.2f}\n", efficiency);
                }
        }
        return 0;
}