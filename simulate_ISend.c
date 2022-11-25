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
 * Student: Jasper Bruin - 12198684
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

    double *cur = malloc((n_local) * sizeof(double));
    double *new = malloc((n_local) * sizeof(double));
    double *old = malloc((n_local) * sizeof(double));

    // malloc array of MPI_Request array to handle the requests with the size of n_local
    MPI_Request *request = malloc(size_Of_Cluster * sizeof(MPI_Request));

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
    } else {
        MPI_Recv(old, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(cur, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }


    for (t = 1; t < t_max; t++) {
        if (left_neightbor != -1) {
            // MPI_Send(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD);
            // MPI_Recv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Isend(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
            MPI_Irecv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
        }

        if (right_neightbor != size_Of_Cluster) {
            //MPI_Send(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD);
            //MPI_Recv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &status);
            MPI_Isend(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
            MPI_Irecv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
        }

        // we calculate the wave equation we are depending on the current and old array we need to wait
        for (int i = 1; i < n_local; i++) {
            MPI_Wait(&request[process_Rank], &status);
            new[i] = 2 * cur[i] - old[i] + 0.15 * (cur[i - 1] - (2 * cur[i] - cur[i + 1]));
        }

        // rotate the buffers
        double *tmp = old;
        old = cur;
        cur = new;
        new = tmp;
    }

    if (process_Rank > 0) {
        for (t = 1; t < size_Of_Cluster; t++) {
            //MPI_Send(cur, n_local, MPI_DOUBLE, 0, t, MPI_COMM_WORLD);
            MPI_Isend(cur, n_local, MPI_DOUBLE, 0, t, MPI_COMM_WORLD, &request[t]);
        }
    } else {
        for (int i = 0; i <= n_local; i++) {
            current_array[i] = cur[i];
        }
        for (int i = 1; i < size_Of_Cluster; i++) {
            //MPI_Recv(current_array + (i * n_local), n_local, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            MPI_Irecv(current_array + (i * n_local), n_local, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request[i]);
        }
    }

    free(cur);
    free(new);
    free(old);
    free(request);
    MPI_Finalize();

    // set first and last element to 0
    current_array[0] = 0;
    current_array[i_max] = 0;

    printf("Returning %d\n", process_Rank);
    return current_array;

}

