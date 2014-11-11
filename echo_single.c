
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
 * SetLandSingle: defines a geometric partitioning of the problem.
 *          and initialises the local rabbit and fox populations.
 * EvolveSingle:  Is called repeatedly to compute the next generation
 *          of both foxes and rabbits. EvolveSingle calls the function
 *          FillBorderSingle to enforce the boundary conditions.
 * GetPopulationSingle: Computes the total population of the specie
 *          it is given as argument.
 * FillBorderSingle: Does the actual halo data swaps for the species it
 *          is given has argument.
 *
 ***********************************************************************/

/* Constants and array bounds */
#include "param.h"

/* A few useful macros. */
#include "macro.h"

/* -- Arrays -- */
/* The arrays storing the number of animals in each stretch of land
 * have borders that are used to store halo data thus simplifying the
 * coding of periodic conditions. The actual data for each animal
 * population is stored in columns 1 through WE_Size and in rows 1 through
 * NS_Size of the arrays below.  from nearest neighbour The halo data is
 * stored in rows 0 and NS_Size+1 and in columns 0 and WE_Size+1 of those
 * arrays.
 */

float Rabbit[NS_Size+2][WE_Size+2];
float Fox[NS_Size+2][WE_Size+2];

/* The next two arrays are used in function EvolveSingle() to compute
 * the next generation of rabbits and foxes.
 */
float TRabbit[NS_Size+2][WE_Size+2];
float TFox[NS_Size+2][WE_Size+2];


/***********************************************************************
 *
 * Initialise the populations of foxes and rabbits.
 *
 ***********************************************************************/
int SetLandSingle ( float Rabbit[NS_Size+2][WE_Size+2],
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
    
    /* Fill the arrays for foxes and rabbits. */
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
 * Set the margin of one of the local data array.
 *
 ***********************************************************************/
int FillBorderSingle(
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
 * Compute the next generation of foxes and rabbits.
 *
 ***********************************************************************/
int EvolveSingle(
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
    err = FillBorderSingle(Rabbit);
    
    /* Fill-in the border of the local Fox array. */
    err = FillBorderSingle(Fox);
    
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
 * Compute the number of individuals in one animal population.
 *
 ***********************************************************************/
int GetPopulationSingle(
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

void echoSingle ()
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
    
    /* Initialise the problem. */
    err = SetLandSingle(Rabbit,Fox,model);
    
    /* Iterate. */
    for( k=1; k<=NITER; k++) {
        err = EvolveSingle(Rabbit,Fox,model);
        if( !(k%PERIOD) ) {
            err = GetPopulationSingle(Rabbit,&nbrab);
            err = GetPopulationSingle(Fox,&nbfox);
            printf("Year %d: %.0f rabbits and %.0f foxes\n",
                   k, nbrab, nbfox);
        }
    }
    
}


