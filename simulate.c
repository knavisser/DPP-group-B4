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

    // print all values:
    if (process_Rank == 0) {
        for (t = 1; t < size_Of_Cluster - 1; t++) {
            for (int i = 1; i < n_local; i++) {
                cur[i] = current_array[i];
            }
            MPI_Send(cur, n_local, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);

            // initialize (old[1:n_local]
            //printf("initialize (old[1:n_local] for process %d \n", t);
            for (int i = 1; i < n_local; i++) {
                old[i] = old_array[i];
            }
            // send (t, old[1:n_nlocal]
            //printf("send (t, old[1:n_nlocal] for process %d \n", t);
            MPI_Send(old, n_local, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);

        }

        // initialize (cur[1:n_local]
        //printf("initialize (cur[1:n_local] for process %d \n", process_Rank);
        for (int i = 1; i < n_local; i++) {
            cur[i] = current_array[i];
        }

        // initialize (old[1:n_local]
        //printf("initialize (old[1:n_local] for process %d \n", process_Rank);
        for (int i = 1; i < n_local; i++) {
            old[i] = old_array[i];
        }

    }
        // Worker processes
    else {
        //printf("receive (t, cur[1:n_nlocal] for process %d \n", process_Rank);
        // cur[1:n_local] = receive ( 0)
        MPI_Recv(cur, n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        //printf("receive (t, old[1:n_nlocal] for process %d \n", process_Rank);
        // old[1:n_local] = receive ( 0)
        MPI_Recv(old, n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    printf("process_Rank: %d, size_Of_Cluster: %d, n_local: %d, left_neightbor: %d, right_neightbor: %d \n", process_Rank, size_Of_Cluster, n_local, left_neightbor, right_neightbor);


//
//    // loop over t-max
//    printf("n_local %d\n", n_local);
//    for (t = 1; t < t_max; t++) {
//        printf("t = %d for process %d\n", t, process_Rank);
//
//        if (left_neightbor != -1) {
//            printf("send (t, cur[1] for process %d\n", process_Rank);
//            fflush(0);
//            MPI_Send(&cur[1], 1, MPI_DOUBLE, left_neightbor, tag, MPI_COMM_WORLD);
//        }
//
//        printf("Receiving cur[n_local] to right_neightbor %d", right_neightbor);
//        MPI_Recv(&cur[n_local + 1], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD, &status);
//
//        printf("Sending cur[n_local] to right_neightbor %d", right_neightbor);
//        if (right_neightbor != size_Of_Cluster) {
//            printf("send (t, cur[n_local] for process %d", process_Rank);
//            MPI_Send(&cur[n_local], 1, MPI_DOUBLE, right_neightbor, tag, MPI_COMM_WORLD);
//        }
//
//        printf("Sending cur[1] to left_neightbor %d", left_neightbor);
//        // cur[0] = receive (left_neighbor)
//        MPI_Recv(&cur[0], 1, MPI_DOUBLE, left_neightbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//        printf("Simulating wave for process %d", process_Rank);
//        for (int i = 1; i < n_local + 1; i++) {
//            new[i] = 2.0 * cur[i] - old[i] + 0.15 * (cur[i - 1] - 2.0 * cur[i] + cur[i + 1]);
//        }
//
//        // rotate the buffers
//        double *temp = old_array;
//        old_array = current_array;
//        current_array = next_array;
//        next_array = temp;
//
//        if (process_Rank > 0) {
//            // send (0, cur[1:n_local])
//            MPI_Send(cur, n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//        }
//        else {
//            // write (file, cur[1:n_local])
//            printf("Process %d: cur[1] = %f, old[1] = %f" , process_Rank, cur[1], old[1]);
//
//            for (int i = 1 ; i < n_local - 1; i++) {
//                //cur[1:n_local] = receive (i)
//                MPI_Recv(cur, n_local, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//                //write (file, cur[1:n_local])
//                printf("Process %d: cur[1] = %f, old[1] = %f" , process_Rank, cur[1], old[1]);
//            }
//        }
//    }

    free(old_array);
    free(current_array);

    // finalize MPI
    MPI_Finalize();

    printf("Finalised");

    return next_array;
}

