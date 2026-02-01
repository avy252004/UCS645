#include <iostream>
#include <vector>
#include <omp.h>
#include <iomanip>

using namespace std;

// Size of the matrix (N x N)
// Start with 500 or 1000. 2000+ might take a while on sequential.
const int N = 1000; 

void initialize(vector<vector<double>>& A) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = (double)rand() / RAND_MAX;
}

int main() {
    vector<vector<double>> A(N, vector<double>(N));
    vector<vector<double>> B(N, vector<double>(N));
    vector<vector<double>> C(N, vector<double>(N, 0.0));

    initialize(A);
    initialize(B);

    cout << "Matrix Size: " << N << "x" << N << endl;
    cout << "Threads: " << omp_get_max_threads() << endl;
    cout << "------------------------------------------" << endl;

    // 1. SEQUENTIAL
    double start = omp_get_wtime();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    double t_seq = omp_get_wtime() - start;
    cout << "Sequential Time:      " << fixed << setprecision(4) << t_seq << "s" << endl;

    // Reset C
    fill(C.begin(), C.end(), vector<double>(N, 0.0));

    // 2. BASIC OPENMP
    start = omp_get_wtime();
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < N; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    double t_par = omp_get_wtime() - start;
    cout << "Basic OpenMP Time:    " << t_par << "s (Speedup: " << t_seq/t_par << "x)" << endl;

    // Reset C
    fill(C.begin(), C.end(), vector<double>(N, 0.0));

    // 3. TRANSPOSED OPENMP (Optimized Memory Access)
    start = omp_get_wtime();
    vector<vector<double>> BT(N, vector<double>(N));
    
    // Transpose B
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            BT[j][i] = B[i][j];
        }
    }

    // Multiply using Transposed B
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < N; k++) {
                sum += A[i][k] * BT[j][k]; // Sequential access for both!
            }
            C[i][j] = sum;
        }
    }
    double t_opt = omp_get_wtime() - start;
    cout << "Optimized (Transp):   " << t_opt << "s (Speedup: " << t_seq/t_opt << "x)" << endl;

    return 0;
}