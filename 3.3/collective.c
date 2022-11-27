// Assignment 3.3: Collective communication
//Collective communication is a form of structured communication, where, instead of a dedicated sender and a dedicated receiver, all MPI processes of a given communicator participate. A simple example is broadcast communication where a given MPI process sends one message to all other MPI processes in the communicator. Write an MPI function in C:
//int MYMPI_Bcast ( void *buffer ,
//                  buffer address
//                  buffer size
//                  datatype of entry root process (sender) commuicator
//
// The function should implement broadcast communication by means of point-to-point communication.
// Each MPI process of the given communicator is assumed to execute the broadcast function with the same message
// parameters and root process argument. The root process is supposed to send the content of its own buffer to all other
// processes of the communicator. Any non-root process uses the given buffer for receiving the corresponding data from the
// root process.A particular advantage of collective communication operations is that their implementation can take
// advantage of the underlying physical network topology without making an MPI-based application program specific to any
// such topology. For your implementation of broadcast assume a 1-dimensional ring topology, where each node only has
// communication links with its two direct neighbours with (circularly) increasing and decreasing MPI process ids.
// While any MPI process may, nonetheless, send messages to any other MPI process, messages need to be routed through a
// number of intermediate nodes. Communication cost can be assumed to be linear in the number of nodes involved. Aim for an efficient implementation of your function on an (imaginary) ring network topology.
// No performance-analysis experiments are required for this assignment. However, you should empirically check that your
// solution actually works, and describe your solution in detail in your report: how it is adapted to the given network
//  topology and how many atomic communication events (sending a message from one node to a neighbouring node) are required
//  to complete one broadcast.
//
// Source: https://mpitutorial.com/tutorials/mpi-broadcast-and-collective-communication/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX 100

void MYMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
{
    // This function differs from the MPI_Send and MPI_Recv functions in that it is collective, meaning that all processes
    // in the communicator must call it. The root process is the process that sends the message to all other processes.
    // Furthermore, the root process specifies the message to be sent, and the other processes receive the message.
    int rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    // the root process sends the message to all other processes
    if (rank == root){
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

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
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
