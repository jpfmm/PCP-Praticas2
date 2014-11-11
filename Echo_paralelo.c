
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
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "echo_single.h"
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

/*** Tamanhos dos offSets em variaveis globais  ***/
int offsetNS, offsetWE, coords[2], rank, k;



/***********************************************************************
 *
 * Initialise the populations of foxes and rabbits.
 *
 ***********************************************************************/
int SetLand ( float Rabbit[offsetNS+2][offsetWE+2], float Fox[offsetNS+2][offsetWE+2], float model[2][3], int landNS, int landWE)
{
    int err;
    int gi, gj, i, j;
    
    err = 0;
    
    /* Set the parameters of the predator-prey model. */
    model[RABBIT][SAME] = -0.2;
    model[RABBIT][OTHER] = 0.6;
    model[RABBIT][MIGRANT] = 0.01;
    model[FOX][OTHER] = 0.6;
    model[FOX][SAME] = -1.8;
    model[FOX][MIGRANT] = 0.02;
    
    /* Fill the arrays for foxes and rabbits. */
    gj=landWE*coords[1]+1;
    for( j=1; j<=offsetWE; j++) {
        gi=landNS*coords[0]+1;
        for( i=1;  i<=offsetNS; i++) {
            Rabbit[i][j] =
            128.0*(gi-1)*(NS_Size-gi)*(gj-1)*(WE_Size-gj) /
            (float)(NS_Size*NS_Size*WE_Size*WE_Size);
            Fox[i][j] =
            8.0*(gi/(float)(NS_Size)-0.5)*(gi/(float)(NS_Size)-0.5)+
            8.0*(gj/(float)(WE_Size)-0.5)*(gj/(float)(WE_Size)-0.5);
            gi++;
        }
        gj++;
    }
    
    /* Return the error code. */
    return(err);
    
}

/***********************************************************************
 * 
 * Compute the next generation of foxes and rabbits.
 * 
 ***********************************************************************/
int Evolve(
float        Rabbit[offsetNS+2][offsetWE+2],
float        Fox[offsetNS+2][offsetWE+2],
float        TRabbit[offsetNS+2][offsetWE+2],
float        TFox[offsetNS+2][offsetWE+2],
float        model[2][3])
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

    
    /*  Não é preciso fill board porque este já é feito na comunicação feita antes de entrar nesta função

    // Update the local population data. */
    for( gj=1; gj<=offsetWE; gj++){
        for( gi=1; gi<=offsetNS; gi++){
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
    for( gj=1; gj<=offsetWE; gj++){
        for( gi=1; gi<=offsetNS; gi++){
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
int GetPopulation(
float    Animal[offsetNS+2][offsetWE+2],
float    *tcount)
{
    int   err;
    int   i, j;
    float p;

    err = 0;
    
    /* Sum population. */
    p = 0.0;
    for( j=1; j<=offsetWE; j++)
        for( i=1; i<=offsetNS; i++){
            p = p + Animal[i][j];
        }

    *tcount = p;

    /* Return the error code. */
    return(err);
}


int main (int argc, char *argv[])
{
    MPI_Request reqSR[4], reqRR[4], reqSF[4], reqRF[4];
    MPI_Status statRR[4], statRF[4], statSR[4], statSF[4];
    MPI_Comm cartcomm;
    
    
    int n_proc, nbrs[4], dims[2], periods[2]={1,1}, reorder=1;
    int landNS, landWE, err,i;
    float sumFox, sumRabb, nbrab, nbfox, model[2][3];
    double time;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    
    if(rank==0){
        time= MPI_Wtime();
        printf("N_proc:%d",n_proc);
    }
    
    
    /****************************************************
     **********    CASO DE 1 PROCESSO  ******************
     ***************************************************/
    if (n_proc==1) {
        echoSingle();
    }else{
        
        /****************************************************
         **********+++    MULTI PROCESSOS  ******************
         ***************************************************/
    
        int lado = sqrt(n_proc);
        dims[0] = lado;
        dims[1] = lado;
        
        if((lado * lado) != n_proc){
            if(rank==0)
                printf("ERRO: Numero incorreto de processos\n");
            MPI_Finalize();
            exit(0);
        }
        
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
        MPI_Comm_rank(cartcomm, &rank);
        MPI_Cart_coords(cartcomm, rank, 2, coords);
        MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
        MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
        
        
        //Actualizar offsets de cada processo
        landNS = offsetNS = NS_Size / lado;
        landWE = offsetWE = WE_Size / lado;
        
        if(coords[0] == (lado-1)){
            offsetNS += NS_Size % lado;
        }
        if(coords[1] == (lado-1)){
            offsetWE += WE_Size % lado;
        }
        
        //Buffers para envio e receção de dados
        float buf_sendFoxN[offsetWE],buf_sendFoxS[offsetWE],buf_sendFoxW[offsetNS],buf_sendFoxE[offsetNS];
        float buf_recvFoxN[offsetWE],buf_recvFoxS[offsetWE],buf_recvFoxW[offsetNS],buf_recvFoxE[offsetNS];
        float buf_sendRabbitN[offsetWE],buf_sendRabbitS[offsetWE],buf_sendRabbitW[offsetNS],buf_sendRabbitE[offsetNS];
        float buf_recvRabbitN[offsetWE],buf_recvRabbitS[offsetWE],buf_recvRabbitW[offsetNS],buf_recvRabbitE[offsetNS];
        
        float Rabbit[offsetNS+2][offsetWE+2];
        float Fox[offsetNS+2][offsetWE+2];
        
        /* The next two arrays are used in function Evolve() to compute
         * the next generation of rabbits and foxes.
         */
        float TRabbit[offsetNS+2][offsetWE+2];
        float TFox[offsetNS+2][offsetWE+2];
        
        //Inicialização das comunicações
        
        //*********  Raposas   **************
        //Enviar
        //Cima e baixo
        MPI_Send_init(&buf_sendFoxN[0], offsetWE, MPI_FLOAT, nbrs[UP], 0, cartcomm, &reqSF[UP]);
        MPI_Send_init(&buf_sendFoxS[0], offsetWE, MPI_FLOAT, nbrs[DOWN], 0, cartcomm, &reqSF[DOWN]);
        
        //Esquerda e direita
        MPI_Send_init(&buf_sendFoxW[0], offsetNS, MPI_FLOAT, nbrs[LEFT], 0, cartcomm, &reqSF[LEFT]);
        MPI_Send_init(&buf_sendFoxE[0], offsetNS, MPI_FLOAT, nbrs[RIGHT], 0, cartcomm, &reqSF[RIGHT]);
        
        //Receber
        //Cima e Baixo
        MPI_Recv_init(&buf_recvFoxS[0], offsetWE, MPI_FLOAT, nbrs[DOWN], 0, cartcomm, &reqRF[DOWN]);
        MPI_Recv_init(&buf_recvFoxN[0], offsetWE, MPI_FLOAT, nbrs[UP], 0, cartcomm, &reqRF[UP]);
        
        //Esquerda e direita
        MPI_Recv_init(&buf_recvFoxE[0], offsetNS, MPI_FLOAT, nbrs[RIGHT], 0, cartcomm, &reqRF[RIGHT]);
        MPI_Recv_init(&buf_recvFoxW[0], offsetNS, MPI_FLOAT, nbrs[LEFT], 0, cartcomm, &reqRF[LEFT]);
        
        //*********  Coelhos   ***************
        //Enviar
        //Cima e baixo
        MPI_Send_init(&buf_sendRabbitN[0], offsetWE, MPI_FLOAT, nbrs[UP], 0, cartcomm, &reqSR[UP]);
        MPI_Send_init(&buf_sendRabbitS[0], offsetWE, MPI_FLOAT, nbrs[DOWN], 0, cartcomm, &reqSR[DOWN]);
        
        //Esquerda e direita
        MPI_Send_init(&buf_sendRabbitW[0], offsetNS, MPI_FLOAT, nbrs[LEFT], 0, cartcomm, &reqSR[LEFT]);
        MPI_Send_init(&buf_sendRabbitE[0], offsetNS, MPI_FLOAT, nbrs[RIGHT], 0, cartcomm, &reqSR[RIGHT]);
        
        //Receber
        //Cima e Baixo
        MPI_Recv_init(&buf_recvRabbitS[0], offsetWE, MPI_FLOAT, nbrs[DOWN], 0, cartcomm, &reqRR[DOWN]);
        MPI_Recv_init(&buf_recvRabbitN[0], offsetWE, MPI_FLOAT, nbrs[UP], 0, cartcomm, &reqRR[UP]);
        
        //Esquerda e direita
        MPI_Recv_init(&buf_recvRabbitE[0], offsetNS, MPI_FLOAT, nbrs[RIGHT], 0, cartcomm, &reqRR[RIGHT]);
        MPI_Recv_init(&buf_recvRabbitW[0], offsetNS, MPI_FLOAT, nbrs[LEFT], 0, cartcomm, &reqRR [LEFT]);
        
        
        /* Initialise the problem. */
        err = SetLand(Rabbit,Fox,model,landNS, landWE);
        
        // Iterate.
        for( k=1; k<=NITER; k++) {
            
            /******************************************************
             ****    Começa comunicação de actualização    ********
             ******************************************************/
            
            
            //**************  Envios ***************/
            //Raposas
            //Cima e baixo
            for(i=1; i <= offsetWE; i++)
                buf_sendFoxN[i-1] = Fox[1][i];
            MPI_Start(&reqSF[UP]);
            
            for(i=1; i <= offsetWE; i++)
                buf_sendFoxS[i-1] = Fox[offsetNS][i];
            MPI_Start(&reqSF[DOWN]);
            
            //Esquerda e direita
            for(i=1; i <= offsetNS; i++)
                buf_sendFoxW[i-1] = Fox[i][1];
            MPI_Start(&reqSF[LEFT]);
            
            for(i=1; i <= offsetNS; i++)
                buf_sendFoxE[i-1] = Fox[i][offsetWE];
            MPI_Start(&reqSF[RIGHT]);
            
            //Coelhos
            //Cima e baixo
            for(i=1; i <= offsetWE; i++)
                buf_sendRabbitN[i-1] = Rabbit[1][i];
            MPI_Start(&reqSR[UP]);
            
            for(i=1; i <= offsetWE; i++)
                buf_sendRabbitS[i-1] = Rabbit[offsetNS][i];
            MPI_Start(&reqSR[DOWN]);
            
            //Esquerda e direita
            for(i=1; i <= offsetNS; i++)
                buf_sendRabbitW[i-1] = Rabbit[i][1];
            MPI_Start(&reqSR[LEFT]);
            
            for(i=1; i <= offsetNS; i++)
                buf_sendRabbitE[i-1] = Rabbit[i][offsetWE];
            MPI_Start(&reqSR[RIGHT]);
            
            
            //**************  Recepção ***************/
            //Raposas
            //Cima e baixo
            MPI_Start(&reqRF[DOWN]);
            MPI_Start(&reqRF[UP]);
            
            //Esquerda e direita
            MPI_Start(&reqRF[RIGHT]);
            MPI_Start(&reqRF[LEFT]);
            
            //Coelhos
            //Cima e baixo
            MPI_Start(&reqRR[DOWN]);
            MPI_Start(&reqRR[UP]);
            
            //Esquerda e direita
            MPI_Start(&reqRR[RIGHT]);
            MPI_Start(&reqRR[LEFT]);
            
            
            //Esperar pelos Receives e aplicar alterações nos quadros
            //Raposas
            MPI_Waitall(4, reqRR , statRR);
            for(i=1; i <= offsetWE; i++)
                Fox[offsetNS+1][i] = buf_recvFoxS[i-1];
            for(i=1; i <= offsetWE; i++)
                Fox[0][i] = buf_recvFoxN[i-1];
            for(i=1; i <= offsetNS; i++)
                Fox[i][offsetWE+1] = buf_recvFoxE[i-1];
            for(i=1; i <= offsetNS; i++)
                Fox[i][0] = buf_recvFoxW[i-1];
            
            //Coelhos
            MPI_Waitall(4, reqRF, statRF);
            for(i=1; i <= offsetWE; i++)
                Rabbit[offsetNS+1][i] = buf_recvRabbitS[i-1];
            for(i=1; i <= offsetWE; i++)
                Rabbit[0][i] = buf_recvRabbitN[i-1];
            for(i=1; i <= offsetNS; i++)
                Rabbit[i][offsetWE+1] = buf_recvRabbitE[i-1];
            for(i=1; i <= offsetNS; i++)
                Rabbit[i][0] = buf_recvRabbitW[i-1];
            
            
            /******************************************************
             ****    Termina comunicação de actualização    ********
             ******************************************************/
            
            err = Evolve(Rabbit,Fox,TRabbit,TFox,model);
            if( !(k%PERIOD) ) {
                err = GetPopulation(Rabbit,&nbrab);
                err = GetPopulation(Fox,&nbfox);
                
                MPI_Reduce(&nbrab, &sumRabb, 1, MPI_FLOAT, MPI_SUM, 0, cartcomm);
                MPI_Reduce(&nbfox, &sumFox, 1, MPI_FLOAT, MPI_SUM, 0, cartcomm);
                
                //if(rank==0)
                  //  printf("Year %d: %.0f rabbits and %.0f foxes\n", k, sumRabb, sumFox);
            }
            
            
            //Esperar que os Sends estejam concluidos para ter a certeza que que já podemos mexer nos buffers
            //(Não creio de que 100% obrigatório)
            MPI_Waitall(4, reqSR , statSR);
            MPI_Waitall(4, reqSF , statSF);
        }
        if(rank==0)
            printf("Year %d: %.0f rabbits and %.0f foxes\n", k, sumRabb, sumFox);

    }
    
    if(rank==0)
            printf("Time: %f\n",MPI_Wtime()-time);
    
    MPI_Finalize();
    return 0;
}





