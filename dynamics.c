#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<zoltan.h>
#include "particles.h"
#include "decompose.h"

/* This function computes forces between two given particles */
int force(struct particle_data *p1,struct particle_data *p2) {

  float dist;

  /* Compute distance between particles */
  dist = sqrt(pow(p1->position.x-p2->position.x,2)+
              pow(p1->position.y-p2->position.y,2)+
              pow(p1->position.z-p2->position.z,2));

  /* Only short range interactions */
  if(dist<epsilon && dist>0.01) {
    p1->force.x-=G*(p1->position.x-p2->position.x)/pow(dist,3);
    p1->force.y-=G*(p1->position.y-p2->position.y)/pow(dist,3);
    p1->force.z-=G*(p1->position.z-p2->position.z)/pow(dist,3);
  }

  return 0;

}

/* This function computes interactions between all particles in the simulation */
int compute_forces(){

  int i,j,k;
  float dist;
  /* EXERCISE 4 - START */
  struct vector3d min,max;
  /* EXERCISE 4 - END */
  int numprocs;
  int procs[size];
  struct export_list_data export_list[4*lnp];
  int nexp=0;
 
  /* Compute local interaction */
  /* EXERCISE 5 - START */
  #pragma omp parallel for private(i,j) shared(particles,lnp)
  /* EXERCISE 5 - END */
  for(i=0;i<lnp;i++) {
    for(j=0;j<lnp;j++) {
      if(i==j) continue;
      force(&particles[i],&particles[j]);
    }
  }

  /* EXERCISE 4 - START */
  for(i=0;i<lnp;i++) {
    min.x=particles[i].position.x-epsilon; max.x=particles[i].position.x+epsilon;
    min.y=particles[i].position.y-epsilon; max.y=particles[i].position.y+epsilon;
    min.z=particles[i].position.z-epsilon; max.z=particles[i].position.z+epsilon;

    Zoltan_LB_Box_Assign(ztn,min.x,min.y,min.z,max.x,max.y,max.z,procs,&numprocs);

    for(j=0;j<numprocs;j++) {
      if(procs[j]==rank) continue;
      export_list[nexp].particle=i;
      export_list[nexp].proc=procs[j];
      nexp++;
    }
  }
  /* EXERCISE 4 - END */
  
  compute_remote_forces(nexp,export_list);

  return 0;

}

/* This function moves particles with respect to forces acting on them */
int move() {

  int i;

  for(i=0;i<lnp;i++) {

    particles[i].position.x+=particles[i].force.x*0.01;
    particles[i].position.y+=particles[i].force.y*0.01;
    particles[i].position.z+=particles[i].force.z*0.01;

    /* Simple periodic boundaries */ 
    if(particles[i].position.x<0) particles[i].position.x=boxsize;
    if(particles[i].position.y<0) particles[i].position.y=boxsize;
    if(particles[i].position.z<0) particles[i].position.z=boxsize;

    if(particles[i].position.x>boxsize) particles[i].position.x=0.0;
    if(particles[i].position.y>boxsize) particles[i].position.y=0.0;
    if(particles[i].position.z>boxsize) particles[i].position.z=0.0;

    particles[i].force.x=0.0;
    particles[i].force.y=0.0;
    particles[i].force.z=0.0;

  }

  return 0;

}

/* This function computes interactions between imported and local particles */ 
int compute_remote_forces(int nexp,struct export_list_data *export_list) {

  int i,j,k;
  struct particle_data export[size][2*lnp];
  struct particle_data import[size][2*lnp];
  int nexport[size];
  int nimport[size];
  float dist;
  MPI_Status status;

  if(nexp==0) return 0;

  for(i=0;i<size;i++) {
    nexport[i]=0;
    nimport[i]=0;
  }

  /* Prepare export buffers */
  for(i=0;i<nexp;i++) {
    int p=export_list[i].proc;
    export[p][nexport[p]]=particles[export_list[i].particle];
    export[p][nexport[p]].force.x=0.0;
    export[p][nexport[p]].force.y=0.0;
    export[p][nexport[p]].force.z=0.0;
    nexport[p]++;
  }

  /* Send information about message sizes */
  MPI_Alltoall(nexport,1,MPI_INT,nimport,1,MPI_INT,MPI_COMM_WORLD);

  /* Send and receive particles */
  for(i=0;i<size;i++) {
    if(i==rank) continue;
    MPI_Sendrecv(&export[i],nexport[i]*sizeof(struct particle_data),MPI_BYTE,i,rank,&import[i],nimport[i]*sizeof(struct particle_data),MPI_BYTE,i,i,MPI_COMM_WORLD,&status);
  }

  /* Compute remote interaction */
  for(i=0;i<size;i++) {
    /* EXERCISE 5 - START */
    #pragma omp parallel for private(j,k) shared(nimport,import,particles,lnp,i)
    /* EXERCISE 5 - END */
    for(j=0;j<nimport[i];j++) {
      for(k=0;k<lnp;k++) {
	force(&import[i][j],&particles[k]);
      }
    }
  }
 
  /* Send and receive results */
  for(i=0;i<size;i++) {
    if(i==rank) continue;
    MPI_Sendrecv(&import[i],nimport[i]*sizeof(struct particle_data),MPI_BYTE,i,rank,&export[i],nexport[i]*sizeof(struct particle_data),MPI_BYTE,i,i,MPI_COMM_WORLD,&status);
  }
 
  for(i=0;i<size;i++) {
    nexport[i]=0;
  }

  /* Update forces */
  for(i=0;i<nexp;i++) {
    int p=export_list[i].proc;
    particles[export_list[i].particle].force.x+=export[p][nexport[p]].force.x;
    particles[export_list[i].particle].force.y+=export[p][nexport[p]].force.y;
    particles[export_list[i].particle].force.z+=export[p][nexport[p]].force.z;
    nexport[p]++;
  }

  return 0;

}
