#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ARRAY_SIZE 10000000 

void MyBcast(double* buffer, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        for (int i = 0; i < size; i++) {
            if (i != root) {
                MPI_Send(buffer, count, MPI_DOUBLE, i, 0, comm);
            }
        }
    } else {
        MPI_Recv(buffer, count, MPI_DOUBLE, root, 0, comm, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double* buffer = (double*)malloc(ARRAY_SIZE * sizeof(double));
    if (rank == 0) {
        for (int i = 0; i < ARRAY_SIZE; i++) buffer[i] = 1.0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_my = MPI_Wtime();
    MyBcast(buffer, ARRAY_SIZE, 0, MPI_COMM_WORLD);
    double end_my = MPI_Wtime();

    if (rank != 0) {
        for (int i = 0; i < ARRAY_SIZE; i++) buffer[i] = 0.0;
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    double start_mpi = MPI_Wtime();
    MPI_Bcast(buffer, ARRAY_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double end_mpi = MPI_Wtime();

    if (rank == 0) {
        double my_time = end_my - start_my;
        double mpi_time = end_mpi - start_mpi;
        printf("--- Broadcast (P=%d) ---\n", size);
        printf("MyBcast Time: %f seconds\n", my_time);
        printf("MPI_Bcast Time: %f seconds\n", mpi_time);
        printf("Performance Gain: %.2fx\n\n", my_time / mpi_time);
    }

    free(buffer);
    MPI_Finalize();
    return 0;
}