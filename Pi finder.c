/*Software Development in the UNIX Environment
Sample C Program

Example C Program to Compute PI Using A Monte Carlo Method

Source code:*/

/* Program to compute Pi using Monte Carlo methods */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define SEED 35791246

int main(int argc, char**){
   double pi, niter, count_local, count, x, y, z;
   int rank_W, size_W;

   /* initialize random numbers */
   srand(SEED);
   count = 0;
   
   MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_W);
	MPI_Comm_size(MPI_COMM_WORLD, &size_W);
   
   if(rank_W == 0){
	printf("Enter the number of iterations used to estimate pi: ");
	scanf("%lg",&niter);
   
	MPI_Bcast(&niter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Reduce(&count_local, &count, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
	
	pi=count/niter*4;
	printf("# of trials= %d , estimate of pi is %g \n",niter,pi);
   }else{
	MPI_Bcast(&niter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	count_local = 0;
	for ( i=0; i<niter; i++) {
      x = (double)rand()/RAND_MAX;
      y = (double)rand()/RAND_MAX;
      z = x*x+y*y;
      if (z<=1) count_local++;
    }
	 
	MPI_Reduce(&count_local, &count, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
	}
   
}