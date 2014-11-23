/***********************

Conway Game of Life

serial version

************************/

#include <stdio.h>

#include <stdlib.h>

#define NI 200        /* array sizes */

#define NJ 200

#define NSTEPS 500    /* number of time steps */

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads()0
#define omp_get_thread_num()0
#endif

int main(int argc, char *argv[]) {

  int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
  int **old, **new;  
  float x;

  if (argc != 2) printf("Argumentos inválidos"); 
  int thread_count = strtol(argv[1], NULL, 10);
  omp_set_num_threads(thread_count);
  
  /* allocate arrays */

  ni = NI + 2;  /* add 2 for left and right ghost cells */
  nj = NJ + 2;
  old = malloc(ni*sizeof(int*));
  new = malloc(ni*sizeof(int*));

  for(i=0; i<ni; i++){
    old[i] = malloc(nj*sizeof(int));
    new[i] = malloc(nj*sizeof(int));
  }

/*  initialize elements of old to 0 or 1 */

  for(i=1; i<=NI; i++){
    for(j=1; j<=NJ; j++){
      x = rand()/((float)RAND_MAX + 1);
      if(x<0.5){
	old[i][j] = 0;
      } else {
	old[i][j] = 1;
      }
    }
  }

  /*  time steps */
  for(n=0; n<NSTEPS; n++){

    /* corner boundary conditions */
    old[0][0] = old[NI][NJ];
    old[0][NJ+1] = old[NI][1];
    old[NI+1][NJ+1] = old[1][1];
    old[NI+1][0] = old[1][NJ];

    /* left-right boundary conditions */

    for(i=1; i<=NI; i++){
      old[i][0] = old[i][NJ];
      old[i][NJ+1] = old[i][1];
    }

    /* top-bottom boundary conditions */

    for(j=1; j<=NJ; j++){
      old[0][j] = old[NI][j];
      old[NI+1][j] = old[1][j];
    }

	#pragma omp parallel for private(im,ip,jm,jp,nsum,j) shared(new,old)
    for(i=1; i<=NI; i++){
      for(j=1; j<=NJ; j++){
	im = i-1;
	ip = i+1;
	jm = j-1;
	jp = j+1;

	nsum =  old[im][jp] + old[i][jp] + old[ip][jp]
	  + old[im][j ]              + old[ip][j ] 
	  + old[im][jm] + old[i][jm] + old[ip][jm];

	switch(nsum){

	case 3:
	  new[i][j] = 1;
	  break;

	case 2:
	  new[i][j] = old[i][j];
	  break;

	default:
	  new[i][j] = 0;
	}
      }
    }
	#pragma omp barrier
    /* copy new state into old state */

    for(i=1; i<=NI; i++){
      for(j=1; j<=NJ; j++){
	old[i][j] = new[i][j];
      }
    }
  }

  /*  Iterations are done; sum the number of live cells */
  isum = 0;
  #pragma omp parallel for private(i) reduction(+:isum)
  for(i=1; i<=NI; i++){
    for(j=1; j<=NJ; j++){
      isum = isum + new[i][j];
    }
  }
  
  #pragma omp master
  printf("\nNumber of live cells = %d\n", isum);

  return 0;
}
