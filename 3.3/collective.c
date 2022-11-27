// Assignment 3.3: Collective communication
// Source: https://mpitutorial.com/tutorials/mpi-broadcast-and-collective-communication/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void MYMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
{
    // This function differs from the MPI_Send and MPI_Recv functions in that it is collective, meaning that all processes
    // in the communicator must call it. The root process is the process that sends the message to all other processes.
    // Furthermore, the root process specifies the message to be sent, and the other processes receive the message.
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    // the root process sends the message to all other processes
    if (rank == root) {
        // send the message to all other processes in the communicator
        for (int i = 0; i < size; i++) {
            // don't send to yourself (root) because you already have the message
            if (i != rank) {
                // send the message to process i in the communicator.
                // i is the rank of the process in the communicator
                MPI_Send(buffer, count, datatype, i, 0, communicator);
            }
        }
    }
    // all other processes receive the message from the root process
    else
    {
        // receive the message from the root process in the communicator
        MPI_Recv(buffer, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
    }
}


int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int data;

    if (world_rank == 0) {
        data = 10000;
        MYMPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Process 0 broadcasting msg %d\n", data);
    } else {
        MYMPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Process %d received msg %d from root process\n", world_rank, data);
    }

    MPI_Finalize();
}

// Output:
//  Process 0 broadcasting msg 10000
//  Process 1 received msg 10000 from root process
//  Process 2 received msg 10000 from root process
//  Process 3 received msg 10000 from root process
