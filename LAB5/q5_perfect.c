#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int is_perfect(int n) {
    if (n <= 1) return 0;
    int sum = 1;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            sum += i;
            if (i * i != n)
                sum += n / i;
        }
    }
    return (sum == n);
}

int main(int argc, char** argv) {
    int rank, size;
    int MAX_VAL = 10000;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("DEBUG: Rank=%d Size=%d\n", rank, size);

    if (rank == 0) {
        double start = MPI_Wtime();

        int next = 2;
        int active = size - 1;
        int perfect_count = 0;
        int response;
        MPI_Status status;

        // 🔥 SEND INITIAL WORK
        for (int i = 1; i < size && next <= MAX_VAL; i++) {
            MPI_Send(&next, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            next++;
        }

        while (active > 0) {
            MPI_Recv(&response, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            int slave = status.MPI_SOURCE;

            if (response == 1)
                perfect_count++;

            if (next <= MAX_VAL) {
                MPI_Send(&next, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
                next++;
            } else {
                int term = -1;
                MPI_Send(&term, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
                active--;
            }
        }

        double end = MPI_Wtime();

        printf("--- Perfect Number Search (P=%d) ---\n", size);
        printf("Perfect Numbers Found: %d\n", perfect_count);
        printf("Total Time: %f seconds\n\n", end - start);
    }

    else {
        int num;

        while (1) {
            MPI_Recv(&num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (num == -1) break;

            int result = is_perfect(num);
            MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}