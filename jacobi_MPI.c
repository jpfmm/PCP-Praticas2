/************************************************************
 * Program to solve a finite difference
 * discretization of the screened Poisson equation:
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * with zero Dirichlet boundary condition using the iterative
 * Jacobi method with overrelaxation.
 *
 * RHS (source) function
 *   f(x,y) = -alpha*(1-x^2)(1-y^2)-2*[(1-x^2)+(1-y^2)]
 *
 * Analytical solution to the PDE
 *   u(x,y) = (1-x^2)(1-y^2)
 *
 * Current Version: Christian Iwainsky, RWTH Aachen University
 * MPI C Version: Christian Terboven, RWTH Aachen University, 2006
 * MPI Fortran Version: Dieter an Mey, RWTH Aachen University, 1999 - 2005
 * Modified: Sanjiv Shah,        Kuck and Associates, Inc. (KAI), 1998
 * Author:   Joseph Robicheaux,  Kuck and Associates, Inc. (KAI), 1998
 *
 * Unless READ_INPUT is defined, a meaningful input dataset is used (CT).
 *
 * Input : n     - grid dimension in x direction
 *         m     - grid dimension in y direction
 *         alpha - constant (always greater than 0.0)
 *         tol   - error tolerance for the iterative solver
 *         relax - Successice Overrelaxation parameter
 *         mits  - maximum iterations for the iterative solver
 *
 * On output
 *       : u(n,m)       - Dependent variable (solution)
 *       : f(n,m,alpha) - Right hand side function
 *
 *************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3

double one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega,
                            int localXcount, int localYcount, int coords[2])
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    for (y = 1; y < (localYcount-1); y++)
    {
        fY = yStart + ((localYcount-2)*coords[0]+(y-1))*deltaY;
        for (x = 1; x < (localXcount-1); x++)
        {
            fX = xStart + ((localXcount-2)*coords[1]+(x-1))*deltaX;
            f = -alpha*(1.0-fX*fX)*(1.0-fY*fY)-2.0*(1.0-fX*fX)-2.0*(1.0-fY*fY);
            updateVal = ((SRC(x-1,y)+SRC(x+1,y))*cx +
                         (SRC(x,y-1)+SRC(x,y+1))*cy +
                         SRC(x,y)*cc - f)/cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }
    return sqrt(error)/((maxXCount)*(maxYCount));
}


/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart,
                     int maxXCount, int maxYCount,
                     double *u,
                     double deltaX, double deltaY,
                     double alpha, int localXcount, int localYcount, int coords[2])
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (localYcount-1); y++)
    {
        fY = yStart + ((localYcount-2)*coords[0]+(y-1))*deltaY;
        for (x = 1; x < (localXcount-1); x++)
        {
            fX = xStart + ((localXcount-2)*coords[1]+(x-1))*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((maxXCount)*(maxYCount));
}


int main(int argc, char **argv)
{
    int n, m, mits;
    double alpha, tol, relax;
    double maxAcceptableError;
    double error, error_local;
    double *u, *u_old, *tmp;
    int allocCount;
    int iterationCount, maxIterationCount;

    MPI_Request reqS[4], reqR[4];
    MPI_Status statS[4], statR[4];
    MPI_Comm cartcomm;

    int period[2]={0,0}, n_proc;

    //Inicializar MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&n_proc);
    int lado = sqrt(n_proc);
    int dims[2] = {lado,lado}, nbrs[4], coords[2];
    int rank;

    if(lado*lado != n_proc){
        printf("Numero de processos invalido\n");
        MPI_Finalize();
        return;
    }


#ifdef READ_INPUT
    printf("Input n,m - grid dimension in x,y direction:\n");
    scanf("%d,%d", &n, &m);
    printf("Input alpha - Helmholts constant:\n");
    scanf("%lf", &alpha);
    printf("Input relax - Successive over-relaxation parameter:\n");
    scanf("%lf", &relax);
    printf("Input tol - error tolerance for iterrative solver:\n");
    scanf("%lf", &tol);
    printf("Input mits - Maximum iterations for solver:\n");
    scanf("%d", &mits);
#else
    n = 200;
    m = 200;
    alpha = 0.8;
    relax = 1.0;
    tol = 1e-7;
    mits = 1000;
#endif
    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);


    //Criar grelha
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,period,1,&cartcomm);
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

    //printf("Hi i'm process %d and my neighbours are:\nLEFT:%d, UP:%d, DOWN:%d, RIGHT: %d\n", rank, nbrs[LEFT], nbrs[UP], nbrs[DOWN], nbrs[RIGHT]);

    int offsetX = n/lado + 2;
    int offsetY = m/lado + 2;

    if(coords[0] == (lado-1)){
        offsetX += n%lado;
    }
    if(coords[1] == (lado-1)){
        offsetY += m%lado;
    }

    //Criar tipo para envio
    //Criar buffers de comunicacao
    double buf_sendXU[offsetX], buf_recvXU[offsetX];
    double buf_sendXD[offsetX], buf_recvXD[offsetX];
    double buf_sendYL[offsetY], buf_recvYL[offsetY];
    double buf_sendYR[offsetY], buf_recvYR[offsetY];

    //Inicializar comunicacao
    //Enviar
    //Cima e baixo
    if(nbrs[UP]!=MPI_PROC_NULL)
    MPI_Send_init(&buf_sendXU[0], offsetX, MPI_DOUBLE, nbrs[UP], 0, cartcomm, &reqS[UP]);
    if(nbrs[DOWN]!=MPI_PROC_NULL)
    MPI_Send_init(&buf_sendXD[0], offsetX, MPI_DOUBLE, nbrs[DOWN], 1, cartcomm, &reqS[DOWN]);
    //Esquerda e direita
    if(nbrs[LEFT]!=MPI_PROC_NULL)
    MPI_Send_init(&buf_sendYL[0], offsetY, MPI_DOUBLE, nbrs[LEFT], 2, cartcomm, &reqS[LEFT]);
    if(nbrs[DOWN]!=MPI_PROC_NULL)
    MPI_Send_init(&buf_sendYR[0], offsetY, MPI_DOUBLE, nbrs[RIGHT], 3, cartcomm, &reqS[RIGHT]);

    //Receber
    //Cime e baixo
    if(nbrs[DOWN]!=MPI_PROC_NULL)
    MPI_Recv_init(&buf_recvXD[0], offsetX, MPI_DOUBLE, nbrs[DOWN], 0, cartcomm, &reqR[DOWN]);
    if(nbrs[UP]!=MPI_PROC_NULL)
    MPI_Recv_init(&buf_recvXU[0], offsetX, MPI_DOUBLE, nbrs[UP], 1, cartcomm, &reqR[UP]);
    //Esquerda e direita
    if(nbrs[RIGHT]!=MPI_PROC_NULL)
    MPI_Recv_init(&buf_recvYR[0], offsetY, MPI_DOUBLE, nbrs[RIGHT], 2, cartcomm, &reqR[RIGHT]);
    if(nbrs[LEFT]!=MPI_PROC_NULL)
    MPI_Recv_init(&buf_recvYL[0], offsetY, MPI_DOUBLE, nbrs[LEFT], 3, cartcomm, &reqR[LEFT]);

    allocCount = offsetY * offsetX;
    // Those two calls also zero the boundary elements
    u = (double*)calloc(sizeof(double), allocCount);
    u_old = (double*)calloc(sizeof(double), allocCount);
    if (u == NULL || u_old == NULL)
    {
        fprintf(stderr, "Not enough memory for two %ix%i matrices\n", offsetY, offsetX);
        exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    int i;

    iterationCount = 0;
    error_local = HUGE_VAL;
    error = HUGE_VAL;

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //Comunica antes de calcular
        //Envia para cima e baixo
        if(nbrs[UP]!=MPI_PROC_NULL){
			for(i=1; i<offsetX-1; i++){
				buf_sendXU[i]=u_old[i+offsetX];
			}
			MPI_Start(&reqS[UP]);
        }
        if(nbrs[DOWN]!=MPI_PROC_NULL){
			for(i=1; i<offsetX-1; i++){
				buf_sendXD[i]=u_old[(offsetY-2)*offsetX+i];
			}
			MPI_Start(&reqS[DOWN]);
        }
        //Envia para esquerda e direita
        if(nbrs[LEFT]!=MPI_PROC_NULL){
			for(i=1; i<offsetY-1; i++){
				buf_sendYL[i]=u_old[offsetX*i+1];
			}
			MPI_Start(&reqS[LEFT]);
        }
        if(nbrs[RIGHT]!=MPI_PROC_NULL){
			for(i=1; i<offsetY-1; i++){
				buf_sendYR[i]=u_old[(offsetX*(i+1))-2];
			}
			MPI_Start(&reqS[RIGHT]);
        }
        //Recebe de cima e de baixo
        if(nbrs[DOWN]!=MPI_PROC_NULL){
			MPI_Start(&reqR[DOWN]);
			for(i=1; i<offsetX-1; i++){
				u_old[(offsetY-1)*offsetX+i]=buf_recvXD[i];
			}
        }
        if(nbrs[UP]!=MPI_PROC_NULL){
			MPI_Start(&reqR[UP]);
			for(i=1; i<offsetX-1; i++){
				u_old[i]=buf_recvXU[i];
			}
        }

        //Recebe da esquerda e direita
        if(nbrs[RIGHT]!=MPI_PROC_NULL){
			MPI_Start(&reqR[RIGHT]);
			for(i=1; i<offsetY-1; i++){
				u_old[(offsetX*(i+1))-1]=buf_recvYR[i];
			}
        }
        if(nbrs[LEFT]!=MPI_PROC_NULL){
			MPI_Start(&reqR[LEFT]);
			for(i=1; i<offsetY-1; i++){
				u_old[offsetX*i]=buf_recvYL[i];
			}
        }

        //Cada um calcula o seu
        printf("Iteration %i\n", iterationCount);
        error_local = one_jacobi_iteration(xLeft, yBottom,
                                     n, m,
                                     u_old, u,
                                     deltaX, deltaY,
                                     alpha, relax, offsetX, offsetY, coords);
        printf("\tError_local %g in coords (%d,%d), rank=%d\n", error_local,coords[0],coords[1], rank);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;

        // AllReduce para ver se continua
        MPI_Allreduce(&error_local, &error, 1, MPI_DOUBLE, MPI_SUM, cartcomm);

        if(rank==0){
           printf("Total error: %g \n", error);
        }
    }
    if(rank==0)
        printf("Residual %g\n",error);

    double absoluteError;
    // u_old holds the solution after the most recent buffers swap
    double absoluteError_local = checkSolution(xLeft, yBottom,
                                         n, m,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha, offsetX, offsetY, coords);
    MPI_Reduce(&absoluteError_local, &absoluteError, 1, MPI_DOUBLE, MPI_SUM, 0, cartcomm);

    if(rank==0)
        printf("The error of the iterative solution is %g\n", absoluteError);

    MPI_Finalize();
    return 0;
}
