
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
 * Purpose:     Program to view a few iteration of the Fox/Rabbit
 *              simulation.
 *
 * Contents:    C source code.
 *
 ************************************************************************/

#include <stdio.h>
#include <math.h>
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
 * Evolve:  Is called repeatedly to compute the next generations
 *          of both foxes and rabbits. Evolve calls the function
 *          FillBorder to enforce the boundary conditions.
 * GetPopulation: Computes the total population of the species
 *          it is given as argument.
 * FillBorder: Does the actual halo data swaps for the species it
 *          is given has argument.
 *
 ***********************************************************************/

/* Constants and array bounds */
#include "param.h"

/* A few useful macros. */
#include "macro.h"

/* -- Arrays -- */
/* The arrays storing the number of animals in each stretch of land
 * have borders that are used as buffers for data from nearest neighbour
 * processors and to simplify the coding of periodic conditions. The
 * data for each animal population is stored in columns 1 through
 * WE_Size and in rows 1 through NS_Size of the arrays below.
 * The halo data is stored in rows 0 and NS_Size+1 and in columns 0 and
 * WE_Size+1 of those arrays.
 */
float Rabbit[NS_Size+2][WE_Size+2];
float Fox[NS_Size+2][WE_Size+2];

/* The next two arrays are used in function Evolve() to compute
 * the next generations of rabbits and foxes.
 */
float TRabbit[NS_Size+2][WE_Size+2];
float TFox[NS_Size+2][WE_Size+2];

void main (int argc, char *argv[])
{

    /* -- Local Arrays -- */
    /* The array named model is used to pass the model parameters
     * across procedures.
     */
    float model[2][3];

    /* Loop indices, bounds and population counters */
    int k;
    float nbrab, nbfox;

    int err;
    char filename[80];

    FILE *fp, *popen();

    fp = popen("gnuplot","w");
    fprintf(fp,"set xrange [1:120]\n");
    fprintf(fp,"set yrange [1:120]\n");
    fprintf(fp,"set zrange [1:6]\n");
    fprintf(fp,"set parametric\n");

    /* Initialise the problem. */
    err = SetLand(Rabbit,Fox,model);
 
    /* Iterate. */
    for( k=1; k<=NITER; k++) {
        err = Evolve(Rabbit,Fox,model); 
        if( !(k%PERIOD) ) {
            err = GetPopulation(Rabbit,&nbrab); 
            err = GetPopulation(Fox,&nbfox); 
            printf("Year %d: %.0f rabbits and %.0f foxes\n",
                k, nbrab, nbfox);
        }
        if(k<6) {
            sprintf(filename,"fox%d",k);
            data(Fox,filename);
            fprintf(fp,"splot '%s'\n",filename);
            fprintf(fp,"!rm %s\n", filename);
        }
    }
    fprintf(fp,"quit");
    pclose(fp);
}
/***********************************************************************
 * 
 * Initialise the populations of foxes and rabbits in the same way
 * as the sequential program does.
 * 
 ***********************************************************************/
int SetLand ( float Rabbit[NS_Size+2][WE_Size+2],
              float Fox[NS_Size+2][WE_Size+2],
              float model[2][3])
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

    /* Fill the local arrays for foxes and rabbits. */
    for( gj=1; gj<=WE_Size; gj++) {
        for( gi=1;  gi<=NS_Size; gi++) {
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
int Evolve(
float        Rabbit[NS_Size+2][WE_Size+2],
float        Fox[NS_Size+2][WE_Size+2],
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

    /* Fill-in the border of the local Rabbit array. */
    err = FillBorder(Rabbit); 

    /* Fill-in the border of the local Fox array. */
    err = FillBorder(Fox); 

    /* Update the local population data. */
    for( gj=1; gj<=WE_Size; gj++){
        for( gi=1; gi<=NS_Size; gi++){
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
    for( gj=1; gj<=WE_Size; gj++){
        for( gi=1; gi<=NS_Size; gi++){
            Rabbit[gi][gj] = MAX( 0.0, TRabbit[gi][gj]);
            Fox[gi][gj] = MAX( 0.0, TFox[gi][gj]);
        }
    }

    /* Return the error code. */
    return(err);
}
/***********************************************************************
 * 
 * Set the margin of one of the local data array. 
 * 
 ***********************************************************************/
int FillBorder(
float        Animal[NS_Size+2][WE_Size+2])
{
    int        err;
    int        i, j;

    err = 0;

    /* Eastward and westward data shifts. */
    for( i=1; i<=NS_Size; i++) {
        Animal[i][0] = Animal[i][WE_Size];
        Animal[i][WE_Size+1] = Animal[i][1];
    }

    /* Northward and southward data shifts. */
    for( j=1; j<=WE_Size; j++) {
        Animal[0][j] = Animal[NS_Size][j];
        Animal[NS_Size+1][j] = Animal[1][j];
    }

    /* Return the error code. */
    return(err);
}
/***********************************************************************
 * 
 * Compute the number of individuals in a population.
 * 
 ***********************************************************************/
int GetPopulation(
float    Animal[NS_Size+2][WE_Size+2],
float    *tcount)
{
    int   err;
    int   i, j;
    float p;

    err = 0;
    
    /* Sum population. */
    p = 0.0;
    for( j=1; j<=WE_Size; j++)
      for( i=1; i<=NS_Size; i++)
            p = p + Animal[i][j];

    *tcount = p;

    /* Return the error code. */
    return(err);
}

data(float Animal[NS_Size+2][WE_Size+2],
     char *filename)
{
    int i, j;
    FILE *fp, *fopen();

    fp = fopen(filename,"w");
    for( j=1;j<=WE_Size;j++ )
        for( i=1;i<=NS_Size;i++)
             fprintf(fp,"%d %d %f\n",j, i, Animal[i][j]);
    fclose(fp);
}



