
/*
 *  A set of functions for handling floating point (double) matrices.
 *
 *  written by Sasa Ilijic, silijic@geof.hr, 2000, 2001, 2002
 *
 */

#define MXFUNS_VERSION "MXFUNS 3C"

#include<stdio.h>

void MxError( char *errs, FILE *logfp, void (*errfoo)(void) );

/* 
 *  This function determines how Mx functions behave on error.
 *
 *  Mx functions of type "double **" * will return NULL in case of an error,
 *  and MxWrite will return EOF.
 * 
 *  errs is a string that will be prepended to the error message.
 *  By default this is "MxCaller".
 * 
 *  logfp is a pointer to a file where error message will be written
 *  in addition to stderr where it always goes.  Give NULL if not
 *  interested in this.
 * 
 *  errfoo is a pointer to a function that is to be called before
 *  returning NULL or EOF (see above).  Give NULL in not interested.
 *  Use this to exit the execution on error if that is your choice.
 * 
 */

double **MxLoad( const char *filename, long *vc, long *vlen );

/*
 *  Loads a matrix from file, vc and/or vlen should be initialized to
 *  expected vector-count and vector-length to reject a non conforming
 *  matrixfile, or set to zero in case the size of matrix is not known
 *  in advance.  If we do not know the size in advance vc and vlen
 *  should be initialized to zero, and if a matrix was successfully
 *  loaded vc and vlen are set to the new values.
 * 
 */

int MxWrite( double **x, long vc, long vlen, const char *filename );

/*
 *  Writes a matrix to a file.
 * 
 */

char *MxFormat( char *fmts );

/*
 *  Makes MxWrite use string fmts as the format string for printing
 *  doubles to te matrixfile.
 *  Empty string restores the default format.
 *  The new setting is returned.
 *
 */

double **MxAlloc( long vc, long vlen );

/*
 *  Allocates memory for a matrix.
 * 
 */

void MxFree( double **x, long vc, long vlen );

/*
 *  De-allocates memory used for a matrix.
 * 
 */

double **MxTrnsp( double **x, long vc, long vlen );

/*
 *  Transposes a matrix,
 *  vc and vlen correspond to matrix before transposition.
 * 
 */

double **MxMply( double **x, double **y, long vc1, long vlen1, long vlen2 );

double **MxFlttn( double **x, long vc, long vlen );

/*
 *  Flattens a matrix into asingle vector matrix,
 *  vc and vlen correspond to matrix before flattening.
 * 
 */

void MxCopy( double **from, double **to, long vc, long vlen );

/*
 *  Copies from to.
 *
 */
