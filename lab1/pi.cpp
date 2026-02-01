#include <iostream>
#include <chrono>
#include <omp.h>
#include <vector>

using namespace std;

int main() {
    const long long num_steps = 1e9;
    const double step = 1.0 / num_steps;

    int max_threads = omp_get_max_threads();

    cout << "Max available threads: " << max_threads << "\n\n";
    cout << "Threads\tTime(ms)\tSpeedup\n";

    double serial_time = 0.0;

    for (int threads = 1; threads <= max_threads; threads *= 2) {
        omp_set_num_threads(threads);

        double sum = 0.0;
        auto start = chrono::steady_clock::now();

        #pragma omp parallel for reduction(+:sum) schedule(static)
        for (long long i = 0; i < num_steps; ++i) {
            double x = (i + 0.5) * step;
            sum += 4.0 / (1.0 + x * x);
        }

        double pi = step * sum;

        auto end = chrono::steady_clock::now();
        chrono::duration<double, milli> elapsed = end - start;

        if (threads == 1) {
            serial_time = elapsed.count();
        }

        double speedup = serial_time / elapsed.count();

        cout << threads << "\t"
             << elapsed.count() << "\t"
             << speedup << "\n";
    }

    return 0;
}
