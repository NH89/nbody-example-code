#include <stdio.h>
#include <stdlib.h>
#include <particles.h>
#include <mpi.h>
/* EXERCISE 1 - START */
#include <sprng.h>
/* EXERCISE 1 - END */

#define SEED 985456376

/* This function generates particle data and sets up few simulation parameters */
int generateParticles(int np) {

  int p;
  /* EXERICES 1 - START */
  int gtype; 
  int *stream;
  /* EXERCISE 1 - END */

  /* EXERCISE 1 - START */
  /* Random Number Generator type: DEFAULT_RNG_TYPE, SPRNG_LFG, SPRNG_LCG, SPRNG_LCG64, SPRNG_CMRG, SPRNG_MLFG, SPRNG_PMLCG */
  gtype = DEFAULT_RNG_TYPE;
  stream = init_sprng(gtype,rank,size,SEED,SPRNG_DEFAULT); 
  /* EXERCISE 1 - END */

  boxsize=64.0;
  epsilon=8.0;

  lnp=np/size;

  particles = (struct particle_data*) malloc(lnp*2*sizeof(struct particle_data));

  for(p=0;p<lnp;p++) {

    particles[p].gid=(ZOLTAN_ID_TYPE)(rank*lnp+p);

    /* EXERICES 1 - START */
    particles[p].position.x = sprng(stream)*boxsize;
    particles[p].position.y = sprng(stream)*boxsize;
    particles[p].position.z = sprng(stream)*boxsize;
    /* EXERCISE 1 - END */

    particles[p].velocity.x = 0.0;
    particles[p].velocity.y = 0.0;
    particles[p].velocity.z = 0.0;    

    particles[p].force.x=0.0;
    particles[p].force.y=0.0;
    particles[p].force.z=0.0;

  }

  return 0;

}


