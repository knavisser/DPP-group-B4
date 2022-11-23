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
    MPI_Status status;

    // Tag: Number to distinguish different categories of messages.
    // It is 0 by default because it is not used in this program.
    int tag = 0;

    int rc = MPI_Init(NULL, NULL);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // N tasks
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);

    // This is the task ID.
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    // n_local = npoints / numtasks: n_local represents the number of points
    int n_local = (int) ceil((double) i_max / size_Of_Cluster);

    // double cur[n_local + 2] double new[n_local + 2]  double old[n_local + 2];
    double *cur = malloc((n_local + 2) * sizeof(double));
    double *new = malloc((n_local + 2) * sizeof(double));
    double *old = malloc((n_local + 2) * sizeof(double));

    // special left treatment: left_neighbor = task_id - 1
    // special right treatment: right_neighbor = task_id + 1
    int left_neightbor = process_Rank - 1;
    int right_neightbor = process_Rank + 1;

    // If the left neighbor is out of bounds, set it to MPI_PROC_NULL.
    if (left_neightbor < 0) {
        left_neightbor = MPI_PROC_NULL;
    }

    // If the right neighbor is out of bounds, set it to MPI_PROC_NULL.
    if (right_neightbor > size_Of_Cluster - 1) {
        right_neightbor = MPI_PROC_NULL;
    }

    // Initialize the arrays from 1 to local size
    for (int i = 1; i < n_local + 1; i++) {
        cur[i] = current_array[i - 1];
        new[i] = next_array[i - 1];
        old[i] = old_array[i - 1];
    }

    // Start the simulation.
    for (t = 0; t < t_max; t++) {
        // Send the left halo cell to the left neighbor.
        MPI_Send(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD);

        // Send the right halo cell to the right neighbor.
        MPI_Send(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD);

        // Receive the right halo cell from the right neighbor.
        MPI_Recv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &status);

        // Receive the left halo cell from the left neighbor.
        MPI_Recv(&cur[0], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD, &status);

        // Calculate the new values.
        for (int i = 1; i < n_local + 1; i++) {
            new[i] = 2.0 * cur[i] - old[i] + 0.15 * (cur[i - 1] - 2.0 * cur[i] + cur[i + 1]);
        }

        // Swap the arrays.
        double *temp = old;
        old = cur;
        cur = new;
        new = temp;
    }

    // Copy the values back to the arrays.
    for (int i = 1; i < n_local + 1; i++) {
        current_array[i - 1] = cur[i];
        next_array[i - 1] = new[i];
        old_array[i - 1] = old[i];
    }

    // set first and last element to 0 from next_array
    next_array[0] = 0;
    next_array[i_max - 1] = 0;

    // Free the memory.
    free(cur);
    free(new);
    free(old);


    printf("Finished process!\n");
    // finalize MPI
    MPI_Finalize();

    return next_array;
}

