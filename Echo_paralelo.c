
/************************************************************************
 * This file has been written for the purpose of courses given at the
 * Edinburgh Parallel Computing Centre. It is made freely available with
 * the understanding that every copy of this file must include this
 * header and that EPCC takes no responsibility for the use of the
 * enclosed teaching material.
 *
 * Author:      Joel Malard
 *
 * Contact:     epcc-tec@epcc.ed.ac.uk
 *
 * Purpose:     Bare sequential program for the case study.
 *
 * Contents:    C source code.
 *
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <mpi.h>
/***********************************************************************
 * 
 * This file contains C code for modelling a simple predator-prey model.
 * The two animal populations modelled are rabbits and foxes that live on
 * some piece of land. The animal populations are represented by two two-
 * dimensional arrays of size NS_Size by WE_Size. Precisely, the number of
 * foxes that live in the i,j-th stretch of land (with 1<=i<=NS_Size and
 * 1<=j<=WE_Size) is stored as Fox[i][j]. The corresponding figure for
 * rabbits is stored in Rabbit[i][j]. Boundary conditions are enforced so
 * that the land is actually a torus and can be thought to extend
 * periodically beyond the given bounds NS_Size and WE_Size.  The
 * populations of each stretch is updated at each generation according to
 * some simple rules and the total populations are summed at regular
 * intervals of time. The correspondence between array indices, process
 * coordinates and cardinal directions are shown in the diagrams below
 * with arrows pointing in the direction of increasing value.
 * 
 *
 *       +-----J----->    
 *       |                  
 *       |                      +--> East
 *       I  A[I][J]             |        
 *       |                      V 
 *       |                    South
 *       V                   
 *
 *
 * Procedures and their roles:
 *****************************
 *
 * SetLand: defines a geometric partitioning of the problem.
 *          and initialises the local rabbit and fox populations.
 * Evolve:  Is called repeatedly to compute the next generation
 *          of both foxes and rabbits. Evolve calls the function
 *          FillBorder to enforce the boundary conditions.
 * GetPopulation: Computes the total population of the specie
 *          it is given as argument.
 * FillBorder: Does the actual halo data swaps for the species it
 *          is given has argument.
 *
 ***********************************************************************/

/* Constants and array bounds */
#include "param.h"

/* A few useful macros. */
#include "macro.h"

 #define UP    0
 #define DOWN  1
 #define LEFT  2
 #define RIGHT 3

/* -- Arrays -- */
/* The arrays storing the number of animals in each stretch of land
 * have borders that are used to store halo data thus simplifying the
 * coding of periodic conditions. The actual data for each animal
 * population is stored in columns 1 through WE_Size and in rows 1 through
 * NS_Size of the arrays below.  from nearest neighbour The halo data is
 * stored in rows 0 and NS_Size+1 and in columns 0 and WE_Size+1 of those
 * arrays.
 */

void main (int argc, char *argv[])
{
	MPI_Request reqSR[4], reqRR[4], reqSF[4], reqRF[4];
	MPI_Status statSR[4], statRR[4], statSF[4], statRF[4];
	MPI_Comm cartcomm;

	int n_proc, rank, source, dest, nbrs[4], dims[2], periods[2]={1,1}, reorder=1, coords[2], i;
	 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	int lado = sqrt(n_proc);
	dims[0] = lado; 
	dims[1] = lado;

	if((lado * lado) != n_proc){
		printf("ERRO: Numero incorreto de processos\n");
		MPI_Finalize();
		exit;
	}

	//Cria a grelha cartesiana e descobre os vizinhos
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
	MPI_Comm_rank(cartcomm, &rank);
	MPI_Cart_coords(cartcomm, rank, 2, coords);
	MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
	MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

	//Determina o offset de cada processador
	int offsetNS = NS_Size / lado + 2;
	int offsetWE = WE_Size / lado + 2;
	
	if(coords[0] == (lado-1)){
		offsetNS += NS_Size % lado;
	}
	if(coords[1] == (lado-1)){
		offsetWE += WE_Size % lado;
	}
	
	//Cria bufs para comunicação
	int buf_sendNS[offsetNS-2], buf_recvNS[offsetNS-2];
	int buf_sendWE[offsetWE-2], buf_recvWE[offsetWE-2];

	//Cria as grelhas
	float Rabbit[offsetNS][offsetWE];
	float Fox[offsetNS][offsetWE];

	/* The next two arrays are used in function Evolve() to compute
	 * the next generation of rabbits and foxes.
	 */
	float TRabbit[offsetNS][offsetWE];
	float TFox[offsetNS][offsetWE];

	//Inicialização das comunicações
	//Raposas
	
	//Enviar
	//Cima e baixo
	MPI_Send_init(&buf_sendWE[0], offsetWE-2, MPI_FLOAT, nbrs[UP], 0, cartcomm, &reqSF[UP]);
	MPI_Send_init(&buf_sendWE[0], offsetWE-2, MPI_FLOAT, nbrs[DOWN], 1, cartcomm, &reqSF[DOWN]);
	
	//Esquerda e direita
	MPI_Send_init(&buf_sendNS[0], offsetNS-2, MPI_FLOAT, nbrs[LEFT], 2, cartcomm, &reqSF[LEFT]);
	MPI_Send_init(&buf_sendNS[0], offsetNS-2, MPI_FLOAT, nbrs[RIGHT], 3, cartcomm, &reqSF[RIGHT]);
	
	//Receber
	//Cima e Baixo
	MPI_Recv_init(&buf_recvWE[0], offsetWE-2, MPI_FLOAT, nbrs[DOWN], 0, cartcomm, &reqRF[DOWN]);
	MPI_Recv_init(&buf_recvWE[0], offsetWE-2, MPI_FLOAT, nbrs[UP], 1, cartcomm, &reqRF[UP]);
	
	//Esquerda e direita
	MPI_Recv_init(&buf_recvNS[0], offsetNS-2, MPI_FLOAT, nbrs[RIGHT], 2, cartcomm, &reqRF[RIGHT]);
	MPI_Recv_init(&buf_recvNS[0], offsetNS-2, MPI_FLOAT, nbrs[LEFT], 3, cartcomm, &reqRF[LEFT]);
		
	//Coelhos
	
	//Enviar
	//Cima e baixo
	MPI_Send_init(&buf_sendWE[0], offsetWE-2, MPI_FLOAT, nbrs[UP], 4, cartcomm, &reqSR[UP]);
	MPI_Send_init(&buf_sendWE[0], offsetWE-2, MPI_FLOAT, nbrs[DOWN], 5, cartcomm, &reqSR[DOWN]);
	
	//Esquerda e direita
	MPI_Send_init(&buf_sendNS[0], offsetNS-2, MPI_FLOAT, nbrs[LEFT], 6, cartcomm, &reqSR[LEFT]);
	MPI_Send_init(&buf_sendNS[0], offsetNS-2, MPI_FLOAT, nbrs[RIGHT], 7, cartcomm, &reqSR[RIGHT]);
	
	//Receber
	//Cima e Baixo
	MPI_Recv_init(&buf_recvWE[0], offsetWE-2, MPI_FLOAT, nbrs[DOWN], 4, cartcomm, &reqRR[DOWN]);
	MPI_Recv_init(&buf_recvWE[0], offsetWE-2, MPI_FLOAT, nbrs[UP], 5, cartcomm, &reqRR[UP]);
	
	//Esquerda e direita
	MPI_Recv_init(&buf_recvNS[0], offsetNS-2, MPI_FLOAT, nbrs[RIGHT], 6, cartcomm, &reqRR[RIGHT]);
	MPI_Recv_init(&buf_recvNS[0], offsetNS-2, MPI_FLOAT, nbrs[LEFT], 7, cartcomm, &reqRR[LEFT]);
	

	
    /* -- Local Arrays -- */
    /* The array named model is used to pass the model parameters
     * across procedures.
     */
    float model[2][3];

    /* Loop indices, bounds and population counters */
    int k;
    float totFox, totRabb, localFox, localRabb;

    int err;

    /* Initialise the problem. (In each processor)*/
    err = SetLand(Rabbit,Fox,model,offsetNS,offsetWE,coords);
 
    /* Iterate. */
    for( k=1; k<=NITER; k++) {
	
		//Envia e recebe ghosts
		//Raposas
		//Envia
		//Cima e baixo
		for(i=1; i <= offsetWE-2; i++){
			buf_sendWE[i] = Fox[1][i];
		}
		MPI_Start(&reqSF[UP]);
		for(i=1; i <= offsetWE-2; i++){
			buf_sendWE[i] = Fox[offsetNS-2][i];
		}
		MPI_Start(&reqSF[DOWN]);
		//Esquerda e direita
		for(i=1; i <= offsetNS-2; i++){
			buf_sendNS[i] = Fox[i][1];
		}
		MPI_Start(&reqSF[LEFT]);
		for(i=1; i <= offsetWE-2; i++){
			buf_sendNS[i] = Fox[i][offsetWE-2];
		}
		MPI_Start(&reqSF[RIGHT]);
		
		//Recebe
		//Cima e baixo
		MPI_Start(&reqRF[DOWN]);
		for(i=1; i <= offsetWE-2; i++){
			Fox[offsetNS-2][i] = buf_sendWE[i];
		}
		MPI_Start(&reqRF[UP]);
		for(i=1; i <= offsetWE-2; i++){
			Fox[0][i] = buf_sendWE[i];
		}
		//Esquerda e direita
		MPI_Start(&reqRF[RIGHT]);
		for(i=1; i <= offsetNS-2; i++){
			Fox[offsetWE-2][i] = buf_sendNS[i];
		}
		MPI_Start(&reqRF[LEFT]);
		for(i=1; i <= offsetNS-2; i++){
			Fox[0][i] = buf_sendNS[i];
		}
			
		//Coelhos
		//Envia
		//Cima e baixo
		for(i=1; i <= offsetWE-2; i++){
			buf_sendWE[i] = Rabbit[1][i];
		}
		MPI_Start(&reqSR[UP]);
		for(i=1; i <= offsetWE-2; i++){
			buf_sendWE[i] = Rabbit[offsetNS-2][i];
		}
		MPI_Start(&reqSR[DOWN]);
		//Esquerda e direita
		for(i=1; i <= offsetNS-2; i++){
			buf_sendNS[i] = Rabbit[i][1];
		}
		MPI_Start(&reqSR[LEFT]);
		for(i=1; i <= offsetWE-2; i++){
			buf_sendNS[i] = Rabbit[i][offsetWE-2];
		}
		MPI_Start(&reqSR[RIGHT]);
		
		//Recebe
		//Cima e baixo
		MPI_Start(&reqRR[DOWN]);
		for(i=1; i <= offsetWE-2; i++){
			Rabbit[offsetNS-2][i] = buf_sendWE[i];
		}
		MPI_Start(&reqRR[UP]);
		for(i=1; i <= offsetWE-2; i++){
			Rabbit[0][i] = buf_sendWE[i];
		}
		//Esquerda e direita
		MPI_Start(&reqRR[RIGHT]);
		for(i=1; i <= offsetNS-2; i++){
			Rabbit[offsetWE-2][i] = buf_sendNS[i];
		}
		MPI_Start(&reqRR[LEFT]);
		for(i=1; i <= offsetNS-2; i++){
			Rabbit[0][i] = buf_sendNS[i];
		}
		
        err = Evolve(Rabbit,Fox,model,offsetNS,offsetWE,coords); 
        if( !(k%PERIOD) ) {
            err = GetPopulation(Rabbit,&localRabb,offsetNS,offsetWE); 
            err = GetPopulation(Fox,&localFox,offsetNS,offsetWE); 
            
			MPI_Reduce(&localRabb, &totRabb, 1, MPI_FLOAT, MPI_SUM, 0, cartcomm);
			MPI_Reduce(&localFox, &totFox, 1, MPI_FLOAT, MPI_SUM, 0, cartcomm);
			
			if(rank==0){
				printf("Year %d: %.0f rabbits and %.0f foxes\n", k, totRabb, totFox);
			}
		}
    }

}
/***********************************************************************
 * 
 * Initialise the populations of foxes and rabbits.
 * 
 ***********************************************************************/
int SetLand ( float *Rabbit, float *Fox, float model[2][3], int offsetNS, int offsetWE, int coords[2])
{
    int err;
    int gi, gj;

    err = 0;

    /* Set the parameters of the predator-prey model. */
    model[RABBIT][SAME] = -0.2;
    model[RABBIT][OTHER] = 0.6;
    model[RABBIT][MIGRANT] = 0.01;
    model[FOX][OTHER] = 0.6;
    model[FOX][SAME] = -1.8;
    model[FOX][MIGRANT] = 0.02;

    /* Fill the arrays for foxes and rabbits. */
    for( gj=1; gj<=offsetWE-2; gj++) {
        for( gi=1;  gi<=offsetNS-2; gi++) {
            gi += coords[0]*offsetNS-2;
			gj += coords[1]*offsetWE-2;
			Rabbit[gi][gj] =
                128.0*(gi-1)*(NS_Size-gi)*(gj-1)*(WE_Size-gj) /
                (float)(NS_Size*NS_Size*WE_Size*WE_Size);
            Fox[gi][gj] =
                8.0*(gi/(float)(NS_Size)-0.5)*(gi/(float)(NS_Size)-0.5)+ 
                8.0*(gj/(float)(WE_Size)-0.5)*(gj/(float)(WE_Size)-0.5);
        }
    }

/* Return the error code. */
    return(err);

}
/***********************************************************************
 * 
 * Compute the next generation of foxes and rabbits.
 * 
 ***********************************************************************/
int Evolve(float **Rabbit, float **Fox, float model[2][3], int offsetNS, int offsetWE, int coords[2])
{
    int err;
    int gi, gj;
    float AlR, BtR, MuR, AlF, BtF, MuF;

    err = 0;

    AlR = model[RABBIT][SAME];
    BtR = model[RABBIT][OTHER];
    MuR  = model[RABBIT][MIGRANT];
    BtF = model[FOX][SAME];
    AlF = model[FOX][OTHER];
    MuF  = model[FOX][MIGRANT];

    /* Fill-in the border of the local Rabbit array. */
    //err = FillBorder(Rabbit); 

    /* Fill-in the border of the local Fox array. */
    //err = FillBorder(Fox); 

    /* Update the local population data. */
    for( gj=1; gj<=offsetWE-2; gj++){
        for( gi=1; gi<=offsetNS-2; gi++){
            gi += coords[0]*offsetNS-2;
			gj += coords[1]*offsetWE-2;
			TRabbit[gi][gj] = (1.0+AlR-4.0*MuR)*Rabbit[gi][gj] +
                BtR*Fox[gi][gj] +
                MuR*(Rabbit[gi][gj-1]+Rabbit[gi][gj+1]+
                Rabbit[gi-1][gj]+Rabbit[gi+1][gj]);
            TFox[gi][gj] = AlF*Rabbit[gi][gj] +
                (1.0+BtF-4.0*MuF)*Fox[gi][gj] +
                MuF*(Fox[gi][gj-1]+Fox[gi][gj+1]+
                Fox[gi-1][gj]+Fox[gi+1][gj]);
        }
    }

    /* Ensure the numbers in Rabbit and Fox are non-negative. */
    for( gj=1; gj<=offsetWE-2; gj++){
        for( gi=1; gi<=offsetNS-2; gi++){
            gi += coords[0]*offsetNS-2;
			gj += coords[1]*offsetWE-2;
			Rabbit[gi][gj] = MAX( 0.0, TRabbit[gi][gj]);
            Fox[gi][gj] = MAX( 0.0, TFox[gi][gj]);
        }
    }

    /* Return the error code. */
    return(err);
}

/***********************************************************************
 * 
 * Compute the number of individuals in one animal population.
 * 
 ***********************************************************************/
int GetPopulation(float *Animal, float *tcount, int offsetNS, int offsetWE)
{
    int   err;
    int   i, j;
    float p;

    err = 0;
    
    /* Sum population. */
    p = 0.0;
    for( j=1; j<=offsetWE-2; j++)
      for( i=1; i<=offsetNS-2; i++)
            p = p + Animal[i][j];

    *tcount = p;

    /* Return the error code. */
    return(err);
}






