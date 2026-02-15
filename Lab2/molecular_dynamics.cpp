#include <chrono>
#include <iostream>
#include <ostream>
#include <random>
#include <vector>
#include <format>
#include <omp.h>
#define SQUARE(a) (a*a)
#define TOTAL_PARTICLES 1e4

constexpr double EPSILON{1.0};
constexpr double SIGMA{1.0};
constexpr double SIGMASQ{SQUARE(SIGMA)};
struct vec3_t {
        double x{},y{},z{};

        vec3_t() = default;
        vec3_t(double x, double y, double z): x{x}, y{y}, z{z} {}

        vec3_t operator+(const vec3_t& other) const { return {x + other.x, y + other.y, z + other.z}; }
        vec3_t operator-(const vec3_t& other) const { return {x - other.x, y - other.y, z - other.z}; }
        vec3_t& operator+=(const vec3_t& other) {
                x += other.x;
                y += other.y;
                z += other.z;
                return *this;
        }
        vec3_t& operator-=(const vec3_t& other) {
                x -= other.x;
                y -= other.y;
                z -= other.z;
                return *this;
        }
        vec3_t operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
        friend auto operator<<(std::ostream& os, const vec3_t& obj) -> std::ostream&;
};

auto operator<<(std::ostream& os, const vec3_t& obj) -> std::ostream& {
        os << std::format("({:.3f})i + ({:.3f})j + ({:.3f})k\n", obj.x, obj.y, obj.z);
        return os;
}

struct Particles {
        std::vector<vec3_t> position{};
        std::vector<vec3_t> force{};
};

auto get_random_number() -> double {
        static thread_local std::mt19937 mt{std::random_device{}()};
        static thread_local std::uniform_real_distribution<double> range{-5.0, 5.0};
        return range(mt);
}

auto generate_random_vector() -> vec3_t {
        vec3_t random_vector{get_random_number(), get_random_number(), get_random_number()};
        return random_vector;
}

auto init_particles(Particles& particles) -> void {
        int n{static_cast<int>(particles.position.size())};
        #pragma omp parallel for
        for(int i = 0; i<n; ++i){
                particles.position[i] = generate_random_vector();
                particles.force[i] = {0.0, 0.0, 0.0};
        }
};

auto compute_potential_and_force_sequential(Particles& particles, double& total_energy) -> std::chrono::duration<double, std::milli> {
        total_energy = 0.0;
        int n{static_cast<int>(particles.position.size())};

        std::cout << "Using No Threading\n";
        auto no_start{std::chrono::steady_clock::now()};
        for(int i{0}; i<n;++i) {
                vec3_t current_force{0.0, 0.0, 0.0};
                for(int j{0}; j<n; ++j) {
                        if(i == j) continue;
                        vec3_t delta{particles.position[i] - particles.position[j]};

                        double r2{SQUARE(delta.x) + SQUARE(delta.y) + SQUARE(delta.z)};

                        if (r2 < 1e-10) continue;
                        double r2_inv{1.0/r2};
                        double s2_inv{SIGMASQ*r2_inv};
                        double s6_inv{SQUARE(s2_inv) * s2_inv};
                        double s12_inv{SQUARE(s6_inv)};

                        double pair_energy{4 * EPSILON * (s12_inv - s6_inv)};
                        total_energy += pair_energy;

                        double force_scalar{(24.0 * EPSILON * r2_inv) * (2.0 * s12_inv - s6_inv)};
                        vec3_t force_vec{delta * force_scalar};

                        current_force += force_vec;
                }
                particles.force[i] = current_force;
        }
        auto no_end{std::chrono::steady_clock::now()};
        std::chrono::duration<double, std::milli> ms{no_end - no_start};
        std::cout << std::format("Execution time: {:.2f}ms\n", ms.count());
        std::cout << std::format("Total Energy: {}\n", total_energy * 0.5);
        return ms;
}

auto compute_potential_and_force_parallel(Particles& particles, double& total_energy, int num_threads) -> std::chrono::duration<double, std::milli> {
        total_energy = 0.0;
        int n{static_cast<int>(particles.position.size())};

        #pragma omp parallel for
        for(int i = 0; i<n; ++i){
                particles.force[i] = {0.0, 0.0, 0.0};
        }

        total_energy = 0.0;
        std::cout << std::format("With {} threads\n", num_threads);
        auto th_start{std::chrono::steady_clock::now()};
        #pragma omp parallel for reduction(+:total_energy) schedule(dynamic) num_threads(num_threads)
        for(int i = 0; i<n;++i) {
                vec3_t current_force{0.0, 0.0, 0.0};
                for(int j = 0; j<n; ++j) {
                        if(i == j) continue;
                        vec3_t delta{particles.position[i] - particles.position[j]};

                        double r2{SQUARE(delta.x) + SQUARE(delta.y) + SQUARE(delta.z)};

                        if (r2 < 1e-10) continue;
                        double r2_inv{1.0/r2};
                        double s2_inv{SIGMASQ*r2_inv};
                        double s6_inv{SQUARE(s2_inv) * s2_inv};
                        double s12_inv{SQUARE(s6_inv)};

                        double pair_energy{4 * EPSILON * (s12_inv - s6_inv)};
                        total_energy += pair_energy;

                        double force_scalar{(24.0 * EPSILON * r2_inv) * (2.0 * s12_inv - s6_inv)};
                        vec3_t force_vec{delta * force_scalar};

                        current_force += force_vec;
                }
                particles.force[i] = current_force;
        }
        auto th_end{std::chrono::steady_clock::now()};
        std::chrono::duration<double, std::milli> ms{th_end - th_start};
        std::cout << std::format("Execution time: {:.2f}ms\n", ms.count());
        std::cout << std::format("Total Energy: {}\n", total_energy * 0.5);
        return ms;
}

auto main(void) -> int {
        Particles particles{};
        particles.position.resize(TOTAL_PARTICLES);
        particles.force.resize(TOTAL_PARTICLES);
        init_particles(particles);
        int max_threads{omp_get_max_threads()};
        double total_energy{0.0};
        auto seq_ms{compute_potential_and_force_sequential(particles, total_energy)};
        std::cout << '\n';
        for(int num_threads{2}; num_threads <= max_threads; num_threads += 2) {
                auto par_ms{compute_potential_and_force_parallel(particles, total_energy, num_threads)};
                double speed_up{seq_ms.count()/par_ms.count()};
                std::cout << std::format("Speed Up: {:.2f}x\n", speed_up);
                std::cout << std::format("Throughput: {:.2f}\n", TOTAL_PARTICLES/par_ms.count());
                std::cout << std::format("Efficiency: {:.2f}\n\n", speed_up/num_threads);
        }
        return 0;
}