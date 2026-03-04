/*
 * (C) Copyright 2021- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <utmpx.h>
#include <sched.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>

#ifndef NOMPI
#include <mpi.h>
#else
#warning Not compiled with MPI support
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#warning Not compiled with OpenMP support
int omp_get_max_threads(void) {
  return 1;
}
int omp_get_thread_num(void) {
  return 0;
}
#endif


#if defined(__APPLE__)
  // These variables and functions are Linux specific
  typedef struct cpu_set_t cpu_set_t;
  struct cpu_set_t {};
  static int CPU_ISSET(int index, const cpu_set_t* cpu_set) { return 0; }
  static void CPU_ZERO(cpu_set_t* cpu_set) {}
  static int sched_getaffinity(pid_t pid, size_t size, cpu_set_t* cpu_set) { return 0; }
#endif

#define maxprocs 256
#define maxthreads 256

void printbind( int rank, int nthreads, int counts[maxthreads], int procs[maxprocs][maxthreads], char hostname[100] )
{
  /* Print out the binding for one rank */
  int thread, proc ;
  printf("Rank %4d on %16s has %3d threads on cores: ",rank, hostname, nthreads);
  
  for(thread=0;thread<nthreads; thread++)
  {
    /* For each thread print out a compressed list of cores that it is bound to */
    printf("(%d",procs[0][thread]);
    if(procs[1][thread]==procs[0][thread]+1)
    {
      printf("-");
    }
    for(proc=1; proc<counts[thread]; proc++)
    {
      if(procs[proc-1][thread]<procs[proc][thread]-1)
      {
        printf(",%d",procs[proc][thread]);
        if(procs[proc][thread]!=counts[thread]&&procs[proc+1][thread]==procs[proc][thread]+1)
        {
          printf("-");
        }
      }
      else if(procs[proc-1][thread]==procs[proc][thread]-1)
      {
        if(proc==counts[thread]-1||procs[proc+1][thread]>procs[proc][thread]+1)
        {
          printf("%d",procs[proc][thread]);
        }
      }
      
    }
   printf(")");
  }
  printf("\n");
  
  return;
}

int main( int argc , char **argv )
{
  int rank,size,cpu;
  int i,np,tid,counts[maxthreads];
  cpu_set_t myset;
  int myprocs[maxprocs][maxthreads];
  int nthreads;
  char hostname[100];

#ifndef NOMPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  gethostname(hostname,99);

  nthreads=omp_get_max_threads();
  /* Create list of procs for this rank & thread */
  #pragma omp parallel default(none) private(myset, np, i, tid) shared(rank, myprocs, counts)
  {
    /* Each thread creates a list of cores that it is bound to */
    tid=omp_get_thread_num();
    np=sysconf(_SC_NPROCESSORS_ONLN);

    /* Initialise */
    counts[tid]=0;
    CPU_ZERO(&myset);
    for(i=0; i<maxprocs; i++)
    {
      myprocs[i][tid]=-1;
    }

    sched_getaffinity(0, sizeof(cpu_set_t), &myset );
    for(i=0; i<np; i++)
    {
      if(CPU_ISSET(i, &myset) == 1)
      {
        myprocs[counts[tid]][tid]=i;
        counts[tid]++;
      }
    }
  }

#ifndef NOMPI
  if(rank==0)
  {
    /* Rank 0 prints its own information then receives and prints the binding of every other process in turn */
    MPI_Status status;
    int rprocs[maxprocs][maxthreads];
    int rcounts[maxthreads], rnthreads; 
    char rhostname[100];

    printbind( rank, nthreads, counts, myprocs, hostname );

    for(i=1; i<size; i++)
    {
      MPI_Recv( &rnthreads, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status );
      MPI_Recv( &rcounts, maxthreads, MPI_INT, i, 2, MPI_COMM_WORLD, &status );
      MPI_Recv( &rprocs, maxthreads*maxprocs, MPI_INT, i, 3, MPI_COMM_WORLD, &status);
      MPI_Recv( &rhostname, 100, MPI_CHAR, i, 4, MPI_COMM_WORLD, &status);
      printbind ( i, rnthreads, rcounts, rprocs, rhostname );
    }

    printf("\n");
  }
  else
  {
    /* All other ranks send their binding information to rank 0 */
    MPI_Send( &nthreads, 1, MPI_INT, 0, 1, MPI_COMM_WORLD );
    MPI_Send( &counts, maxthreads, MPI_INT, 0, 2, MPI_COMM_WORLD );
    MPI_Send( &myprocs,maxthreads*maxprocs, MPI_INT, 0, 3, MPI_COMM_WORLD );
    MPI_Send( &hostname, 100, MPI_CHAR, 0, 4, MPI_COMM_WORLD );

  }

  MPI_Finalize();

#else
  printbind( 0, nthreads, counts, myprocs, hostname );
#endif

  return 0;
}
