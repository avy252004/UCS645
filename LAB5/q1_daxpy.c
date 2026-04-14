#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 65536

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a = 2.5;
    double *X = NULL, *Y = NULL, *X_seq = NULL;
    int local_N = N / size;
    double *local_X = (double*)malloc(local_N * sizeof(double));
    double *local_Y = (double*)malloc(local_N * sizeof(double));

    double seq_time = 0.0; 

    if (rank == 0) {
        X = (double*)malloc(N * sizeof(double));
        Y = (double*)malloc(N * sizeof(double));
        X_seq = (double*)malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) {
            X[i] = 1.0;
            Y[i] = 2.0;
            X_seq[i] = X[i];
        }

        double start_seq = MPI_Wtime();
        for (int i = 0; i < N; i++) {
            X_seq[i] = a * X_seq[i] + Y[i];
        }
        double end_seq = MPI_Wtime();

        seq_time = end_seq - start_seq; 
        printf("--- DAXPY (P=%d) ---\n", size);
        printf("Uniprocessor Time: %f seconds\n", seq_time);
    }

    double start_mpi = MPI_Wtime();
    MPI_Scatter(X, local_N, MPI_DOUBLE, local_X, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(Y, local_N, MPI_DOUBLE, local_Y, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double start_comp = MPI_Wtime();
    for (int i = 0; i < local_N; i++) {
        local_X[i] = a * local_X[i] + local_Y[i];
    }
    double end_comp = MPI_Wtime();

    MPI_Gather(local_X, local_N, MPI_DOUBLE, X, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double end_mpi = MPI_Wtime();

    if (rank == 0) {
        double mpi_time = end_mpi - start_mpi;
        double comp_time = end_comp - start_comp;
        double comm_time = mpi_time - comp_time;
        double speedup = seq_time / mpi_time;
        double efficiency = (speedup / size) * 100.0;
        double comm_percent = (comm_time / mpi_time) * 100.0;

        printf("MPI Parallel Time: %f seconds\n", mpi_time);
        printf("Speedup: %.2f\n", speedup); 
        printf("Efficiency: %.2f%%\n", efficiency);
        printf("Comm Overhead: %.2f%%\n\n", comm_percent);

        free(X); free(Y); free(X_seq);
    }

    free(local_X); free(local_Y);
    MPI_Finalize();
    return 0;
}