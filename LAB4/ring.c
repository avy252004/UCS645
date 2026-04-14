#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int value = 100;
    double start_time = MPI_Wtime();

    if (size == 1) {
        value += rank;
        double end_time = MPI_Wtime();
        printf("P=1 Time: %f sec\n", end_time - start_time);
    } else {
        if (rank == 0) {
            printf("Process %d started with value: %d\n", rank, value);
            value += rank; 
            MPI_Send(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&value, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double end_time = MPI_Wtime();
            printf("Process 0 received final wrapped value: %d\n", value);
            printf("Ring P=%d Time: %f sec\n", size, end_time - start_time);
        } else {
            MPI_Recv(&value, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process %d received value: %d\n", rank, value);
            value += rank; 
            int next_rank = (rank + 1) % size;
            MPI_Send(&value, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}