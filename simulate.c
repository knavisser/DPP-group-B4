/*
 * simulate.c
 *
 * assignment 3.1
 *
 * This program simulates a 1-D wave equation using the finite difference method. We are parallelizing this program using
 * MPI which is a message passing interface. This program is written in C and uses the MPI functions.
 * Student: Jasper Bruin - 12198684 & Tony Visser -
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "simulate.h"
#include <math.h>


double *simulate(const int i_max, const int t_max, double *old_array,
                 double *current_array, double *next_array) {

    // Initialize MPI
    int process_Rank, size_Of_Cluster, t;
    int tag = 0;
    MPI_Status status;
    int rc = MPI_Init(NULL, NULL);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // Number of tasks and Task ID
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    // Calculate the number of elements per process
    int n_local = (int) ceil((double) i_max / size_Of_Cluster);

    // Calculate the start and end index for each process, + 2 for the ghost cells.
    double *cur = malloc((n_local + 2) * sizeof(double));
    double *new = malloc((n_local + 2) * sizeof(double));
    double *old = malloc((n_local + 2) * sizeof(double));

    int left_neightbor = process_Rank - 1;
    int right_neightbor = process_Rank + 1;

    // If we are the root process, then we want to initialize the arrays and send them to the other processes.
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
    // If we are not the root process, then we want to receive the arrays from the root process.
    } else {
        MPI_Recv(old, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(cur, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }

    // Calculate the start and end index for each process, + 2 for the ghost cells.
    for (t = 1; t < t_max; t++) {
        // If there is a left neightbor, send the first element to the left neightbor.
        if (left_neightbor != -1) {
            MPI_Send(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD);
            MPI_Recv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // If there is a right neightbor, send the last element to the right neightbor.
        if (right_neightbor != size_Of_Cluster) {
            MPI_Send(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD);
            MPI_Recv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &status);
        }

        // Simulate the wave equation for the current time step.
        for (int i = 1; i <= n_local; i++) {
            new[i] = 2.0 * cur[i] - old[i] + 0.15 * (cur[i - 1] - (2 * cur[i] - cur[i + 1]));
        }

        // rotate the buffers
        double *tmp = old;
        old = cur;
        cur = new;
        new = tmp;

        // Set last element of the array to 0.0
         if (process_Rank == 0) {
            cur[0] = 0;
            old[0] = 0;
            new[0] = 0;
         }

        if (process_Rank == size_Of_Cluster - 1) {
            cur[n_local] = 0;
            old[n_local] = 0;
            new[n_local] = 0;
        }

    }

    // Merge the arrays by sending them to the root process.
    if (process_Rank > 0) {
        MPI_Send(cur, n_local, MPI_DOUBLE, 0, process_Rank, MPI_COMM_WORLD);
    } else {
        // If we are the root process, then we want to receive the arrays from the other processes.
        for (int i = 0; i <= n_local; i++) {
            current_array[i] = cur[i];
        }

        for (int i = 1; i < size_Of_Cluster; i++) {
            MPI_Recv(current_array + (i * n_local), n_local, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        }
    }

    // Finalize MPI and return the current array.
    if (process_Rank == 0) {
        MPI_Finalize();
        free(cur);
        free(new);
        free(old);

        current_array[i_max] = 0.0;
        return current_array;
    }
    // If we are not the root process, then we want to free the memory and return NULL.
    else {
        MPI_Finalize();
        free(cur);
        free(new);
        free(old);
        exit(0);
    }
}


