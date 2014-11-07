
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
 * Purpose:     Some useful macros.
 *
 * Contents:    C include file.
 *
 ************************************************************************/

#define MAX( A, B) (((A)<(B))?(B):(A))
#define MIN( A, B) (((A)>=(B))?(B):(A))
#define ABS( A) (((A)<0.0)?(-(A)):(0))
