#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define TOTAL_ELEMENTS 50000000 
int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double multiplier = 0.0;
    
    if (rank == 0) {
        multiplier = 1.5;
    }
    
    double start_time = MPI_Wtime();
    MPI_Bcast(&multiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double start_comp = MPI_Wtime();
    long local_N = TOTAL_ELEMENTS / size;
    double *local_A = (double*)malloc(local_N * sizeof(double));
    double *local_B = (double*)malloc(local_N * sizeof(double));

    for (long i = 0; i < local_N; i++) {
        local_A[i] = 1.0;
        local_B[i] = 2.0 * multiplier;
    }

    double local_dot = 0.0;
    for (long i = 0; i < local_N; i++) {
        local_dot += local_A[i] * local_B[i];
    }
    double end_comp = MPI_Wtime();

    double final_result = 0.0;
    MPI_Reduce(&local_dot, &final_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        double comp_time = end_comp - start_comp;
        double comm_time = total_time - comp_time;
        
        printf("--- Dot Product (P=%d) ---\n", size);
        printf("Final Dot Product: %f\n", final_result);
        printf("Total Time: %f seconds\n", total_time);
        printf("Comm Overhead: %.2f%%\n\n", (comm_time / total_time) * 100.0);
    }

    free(local_A); free(local_B);
    MPI_Finalize();
    return 0;
}