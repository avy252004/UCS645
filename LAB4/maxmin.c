#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct val_rank { int value; int rank; };

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank * 100);
    int local_max = -1, local_min = 1001;
    
    double start_time = MPI_Wtime();

    if (size > 1) { 
        printf("Process %d generated: ", rank);
        for (int i = 0; i < 10; i++) {
            int num = rand() % 1001; 
            printf("%d ", num);
            if (num > local_max) local_max = num;
            if (num < local_min) local_min = num;
        }
        printf("\n");
    } else {
       
        for (int i = 0; i < 10; i++) {
            int num = rand() % 1001; 
            if (num > local_max) local_max = num;
            if (num < local_min) local_min = num;
        }
    }

    struct val_rank local_max_struct = {local_max, rank};
    struct val_rank global_max_struct;
    struct val_rank local_min_struct = {local_min, rank};
    struct val_rank global_min_struct;

    MPI_Reduce(&local_max_struct, &global_max_struct, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min_struct, &global_min_struct, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        if (size > 1) {
            printf("\n--- Global Results ---\n");
            printf("Global Maximum: %d (Found in Process %d)\n", global_max_struct.value, global_max_struct.rank);
            printf("Global Minimum: %d (Found in Process %d)\n", global_min_struct.value, global_min_struct.rank);
        }
        printf("MaxMin P=%d Time: %f sec\n", size, end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}