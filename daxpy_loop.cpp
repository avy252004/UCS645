#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>

int main() {
    const long N = 1 << 16;
    double a = 2.5;

    std::vector<double> X(N, 1.0), Y(N, 2.0);

    std::ofstream file("daxpy.csv");
    file << "threads,time\n";

    for (int threads = 1; threads <= 16; threads *= 2) {
        omp_set_num_threads(threads);

        double start = omp_get_wtime();

        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            X[i] = a * X[i] + Y[i];
        }

        double end = omp_get_wtime();
        double time = end - start;

        std::cout << "Threads: " << threads
                  << " Time: " << time << "\n";

        file << threads << "," << time << "\n";
    }

    file.close();
    return 0;
}
