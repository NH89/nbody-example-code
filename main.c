#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<particles.h>

int main(int argc,char **argv) {

  int i;
  int iter;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(argc<3) {
    if(rank==0) {
      fprintf(stderr,"Parameters are missing.\n");
    }
    exit(1); 
  } else {
    np = atoi(argv[1]);
    niter = atoi(argv[2]);	
  }

  generateParticles(np);

  decompositionInit(argc,argv);

  for(iter=0;iter<niter;iter++) {

    if(rank==0) printf("Iteration %d\n",iter);

    write_particles(iter);
    decompose();
    compute_forces();
    move();

  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();  

  return 0;
}
