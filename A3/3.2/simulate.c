/*
 * simulate.c
 *
 * assignment 3.2
 *
 * This program simulates a 1-D wave equation using the finite difference method. We are parallelizing this program using
 * MPI which is a message passing interface. This program is written in C and uses the MPI functions.
 *
 * This is the optimised version of 3.1 because we use the MPI_ISend and MPI_IRecv functions. These are non-blocking send
 * and receive functions that allow the program to continue while the message is being sent or received. This is done by
 * using the MPI_Request variable. This variable is used to keep track of the status of the message. The MPI_Wait function
 * is used to wait for the message to be sent or received.
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

    // Calculate the number of elements per process and the start and end index for each process, + 2 for the ghost cells.
    int n_local = (int) ceil((double) i_max / size_Of_Cluster);
    double *cur = malloc((n_local) * sizeof(double));
    double *new = malloc((n_local) * sizeof(double));
    double *old = malloc((n_local) * sizeof(double));

    // Request variables for the MPI_Isend and MPI_Irecv functions
    MPI_Request *request = malloc(size_Of_Cluster * sizeof(MPI_Request));


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
    }
        // If we are not the root process, then we want to receive the arrays from the root process.
    else {
        MPI_Recv(old, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(cur, n_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }


    // Calculate the start and end index for each process.
    for (t = 1; t < t_max; t++) {
        // Calculate the start and end index for each process.
        if (left_neightbor != -1) {
            // Send the first element of the array to the left neightbor.
            MPI_Isend(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
            MPI_Irecv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
        }

        // Calculate the start and end index for each process.
        if (right_neightbor != size_Of_Cluster) {
            MPI_Isend(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
            MPI_Irecv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &request[process_Rank]);
        }

        // Calculate the start and end index for each process and wait for the message to be sent or received.
        for (int i = 1; i < n_local; i++) {
            // Calculate the start and end index for each process. + 2 for the ghost cells.
            // Wait for the message to be sent or received.
            MPI_Wait(&request[process_Rank], &status);
            new[i] = 2 * cur[i] - old[i] + 0.15 * (cur[i - 1] - (2 * cur[i] - cur[i + 1]));
            MPI_Wait(&request[process_Rank], &status);
        }

        // rotate the buffers
        double *tmp = old;
        old = cur;
        cur = new;
        new = tmp;
    }

    // If we are the root process, then we want to receive the arrays from the other processes.
    if (process_Rank > 0) {
        for (t = 1; t < size_Of_Cluster; t++) {
            MPI_Isend(cur, n_local, MPI_DOUBLE, 0, t, MPI_COMM_WORLD, &request[t]);
        }
    }
        // If we are not the root process, then we want to send the arrays to the root process.
    else {
        for (int i = 0; i <= n_local; i++) {
            MPI_Wait(&request[process_Rank], &status);
            current_array[i] = cur[i];
        }
        for (int i = 1; i < size_Of_Cluster; i++) {
            MPI_Irecv(current_array + (i * n_local), n_local, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request[i]);
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
        free(request);

        current_array[i_max] = 0.0;
        return current_array;
    }
        // If we are not the root process, then we want to free the memory and return NULL.
    else {
        MPI_Finalize();
        free(cur);
        free(new);
        free(old);
        free(request);
        exit(0);
    }
}


