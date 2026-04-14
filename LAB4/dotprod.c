#include <mpi.h>
#include <stdio.h>

#define N 8

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int A[N] = {1, 2, 3, 4, 5, 6, 7, 8};
    int B[N] = {8, 7, 6, 5, 4, 3, 2, 1};
    int elements_per_proc = N / size;
    int local_A[8], local_B[8]; 

    double start_time = MPI_Wtime();

    MPI_Scatter(A, elements_per_proc, MPI_INT, local_A, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, elements_per_proc, MPI_INT, local_B, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

    int local_dot = 0;
    for (int i = 0; i < elements_per_proc; i++) {
        local_dot += local_A[i] * local_B[i];
    }

    int global_dot = 0;
    MPI_Reduce(&local_dot, &global_dot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        if(size > 1){
            printf("Parallel Dot Product Result: %d\n", global_dot);
            printf("Expected Result: 120\n");
        }
        printf("DotProd P=%d Time: %f sec\n", size, end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}