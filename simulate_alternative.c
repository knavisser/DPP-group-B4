/*
 * simulate.c
 *
 * assignment 2.1
 *
Reconsider the 1-dimensional wave equation studied in assignments 1.1 and 1.2.
 Write a distributed MPI program in C that uses multiple compute nodes to simulate the wave equation following
 the do- main decomposition approach illustrated below (similar to what was explained during the lecture).
 Each MPI process shall only store the relevant parts of the three amplitude buffers in its local memory.
 Add halo cells as necessary, and exchange their values between time steps as needed. Aim at an efficient
 implementation w.r.t. communication and memory, and use blocking MPI calls to implement the communication
 between processes.
 * Student: Jasper Bruin
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "simulate.h"
#include <math.h>


double *simulate(const int i_max, const int t_max, double *old_array,
                 double *current_array, double *next_array) {

    int process_Rank, size_Of_Cluster, t;
    int tag = 0;
    MPI_Status status;
    int rc = MPI_Init(NULL, NULL);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // N tasks / Task ID
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    int n_local = (int) ceil((double) i_max / size_Of_Cluster);

    double *cur = malloc((n_local + 2) * sizeof(double));
    double *new = malloc((n_local + 2) * sizeof(double));
    double *old = malloc((n_local + 2) * sizeof(double));

    int left_neightbor = process_Rank - 1;
    int right_neightbor = process_Rank + 1;

    if (process_Rank == 0) {
        for (int i = 0; i <= n_local; i++) {
            cur[i] = current_array[i];
            old[i] = old_array[i];
        }

        for (t = 1; t < size_Of_Cluster; t++) {
            for (int i = 0; i <= n_local; i++) {
                cur[i] = current_array[i];
                old[i] = old_array[i];
            }
            MPI_Send(cur, n_local, MPI_DOUBLE, t, tag, MPI_COMM_WORLD);
            MPI_Send(old, n_local, MPI_DOUBLE, t, tag, MPI_COMM_WORLD);
        }
    }

    else {
        MPI_Recv(old, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(cur, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }

    for (t = 1; t < t_max; t++) {
        if (left_neightbor != -1) {
            MPI_Send(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD);
            MPI_Recv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (right_neightbor != size_Of_Cluster) {
            MPI_Send(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD);
            MPI_Recv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &status);
        }

        for (int i = 1; i < n_local + 1; i++) {
            new[i] = 2.0 * cur[i] - old[i] + 0.15 * (cur[i - 1] - 2.0 * cur[i] + cur[i + 1]);
        }

        // rotate the buffers
        double *tmp = old;
        old = cur;
        cur = new;
        new = tmp;
    }

    if (process_Rank > 0) {
        for (int i = 1; i < size_Of_Cluster; i++) {
            for (int j = 0; j < n_local + 1; j++) {
                MPI_Send(&cur[j], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            }
        }
    }

    // master receives all the data from the clusters and puts it in current_array
    if (process_Rank == 0) {
        // master puts its own data in current_array
        for (int i = 0; i < n_local + 1; i++) {
            current_array[i] = cur[i];
        }

        // master receives data from the clusters
        for (int i = 1; i < size_Of_Cluster; i++) {
            for (int j = 0; j < n_local + 1; j++) {
                MPI_Recv(&current_array[j], 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
            }
        }
    }


    MPI_Finalize();
    free(cur);
    free(new);
    free(old);
    return current_array;
}


