
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
 * Purpose:     Internal limits and array bounds.
 *
 * Contents:    C include file.
 *
 ************************************************************************/

/* The number of dimensions of the land. */
#define NADIM 2

/* the dimension of the process mesh. */
#define NPDIM 2

/* the number of iterations */
#define NITER 50

/* the frequency at which the populations are surveyed. */
#define PERIOD 1

/* The following 5 constants are used to index nearest neighbours. */
#define HERE 0
#define EAST 1
#define WEST 2
#define NORTH 3
#define SOUTH 4

/* Global upper bounds for the numbers of land cells in the X and Y dimensions,
 * respectively.
 */
#define NS_Size 120
#define WE_Size 120

/* Upper bounds on the size of the first and second dimension of the local
 * arrays that record animal populations.
 */
#define LocalMaxX 120
#define LocalMaxY 120

/* Labels for the dimensions of various arrays */
#define COLUMN 0
#define ROW 1
#define NS 0
#define WE 1

#define RABBIT 0
#define FOX 1

#define SAME 0
#define OTHER 1
#define MIGRANT 2

/* Other constants */
#define TRUE 1
#define FALSE 0



