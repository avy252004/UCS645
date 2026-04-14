#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ARRAY_SIZE 100

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int elements_per_proc = ARRAY_SIZE / size;
    int *full_array = NULL;
    int *local_array = (int*)malloc(elements_per_proc * sizeof(int));

    if (rank == 0) {
        full_array = (int*)malloc(ARRAY_SIZE * sizeof(int));
        for (int i = 0; i < ARRAY_SIZE; i++) full_array[i] = i + 1;
    }

    double start_time = MPI_Wtime();

    MPI_Scatter(full_array, elements_per_proc, MPI_INT, local_array, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

    int local_sum = 0;
    for (int i = 0; i < elements_per_proc; i++) local_sum += local_array[i];

    int global_sum = 0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("ArraySum P=%d Time: %f sec\n", size, end_time - start_time);
        printf("Global Sum: %d (Expected: 5050)\n", global_sum);
        printf("Average Value: %.2f\n", (float)global_sum / ARRAY_SIZE);
        free(full_array);
    }

    free(local_array);
    MPI_Finalize();
    return 0;
}