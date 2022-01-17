
/*
 *  module mxfuns by Sasa Ilijic, silijic@geof.hr, 2000, 2001, 2002
 *
 *  see mxfuns.h
 *
 *  last revision Feb 2002, GF, Zagreb
 *
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include"mxfuns.h"

static char *mxerrs="MxCaller", *mxerr1="MxFunction", mxerr2[256];
static FILE *mxlogfp = NULL;
static void (*mxerrfoo)(void);

#define PrintError { if ( NULL != mxlogfp ) \
        fprintf( mxlogfp, "\n%s : %s : %s\n", mxerrs, mxerr1, mxerr2 );\
        fprintf( stderr,  "\n%s : %s : %s\n", mxerrs, mxerr1, mxerr2 );}

#define TreatError { if( NULL != mxerrfoo ) mxerrfoo(); else return NULL; }

#define MX_DEFAULT_FORMAT "% E   "
static char *mxfmts=MX_DEFAULT_FORMAT;


void MxError( char *errs, FILE *logfp, void (*errfoo)(void) ) {

    mxerrs = errs;
    mxlogfp = logfp;
    mxerrfoo = errfoo;
}


char *MxFormat( char *fmts ) {

    return  mxfmts = strlen( fmts ) ? fmts : MX_DEFAULT_FORMAT;
}


double **MxLoad( const char *filename, long *vc, long *vlen ) {

    int c, ok;
    long i, j, invc, invlen;
    double **mx;
    char *mxerr0, stringbuff[100];
    FILE *fp;

    mxerr0 = mxerr1;
    mxerr1 = "MxLoad";

    invc = *vc;
    invlen = *vlen;

    if( strcmp( "stdin", filename ) ) {
        if( NULL == ( fp = fopen( filename, "r" ) ) ) {
            sprintf( mxerr2, "failed open to read \"%s\"", filename );
            PrintError; TreatError;
        }
    } else {
        fp = stdin;
    }

    if( 2 != fscanf( fp, "# %ld X %ld", vc, vlen ) ) {
        sprintf( mxerr2, "failed parsing header of \"%s\"", filename);
        PrintError; TreatError;
    }

    if( ( 1 > *vc ) || ( 1 > *vlen ) ) {
        sprintf( mxerr2, "funny header in \"%s\"", filename);
        PrintError; TreatError;
    }

    if( ((0<invc)&&(invc!=*vc)) || ((0<invlen)&&(invlen!=*vlen)) ) {
        sprintf( mxerr2,
            "size of \"%s\" is # %ld X %ld (expected is # %ld X %ld)",
            filename, *vc, *vlen, invc, invlen );
        PrintError; TreatError;
    }

    if( NULL == ( mx = MxAlloc( *vc, *vlen ) ) ) {
        TreatError;
    }

    for ( i = 0 ; i < *vlen ; i++ ) {
        for ( j = 0 ; j < *vc ; j++ ) {
            ok = 0;
            do {
                if ( 1 != fscanf( fp, "%99s", stringbuff ) ) {
                    sprintf( mxerr2,
                        "failed parsing body of \"%s\"", filename );
                    MxFree( mx, *vc, *vlen );
                    PrintError; TreatError;
                } else if ( '#' == *stringbuff ) {
                    do {
                        c = fgetc( fp );
                    } while ( ( EOF != c ) && ( '\n' != c ) );
                } else {
                    if ( 1 != sscanf( stringbuff, "%lg", *(mx+j)+i ) ) {
                        sprintf( mxerr2,
                            "failed parsing body of \"%s\"", filename );
                        MxFree( mx, *vc, *vlen );
                        PrintError; TreatError;
                    } else {
                        ok = 1;
                    }
                }
            } while ( !ok );
        }
    }

    fclose(fp);

    mxerr1 = mxerr0;

    return mx;
}


int MxWrite( double **x, long vc, long vlen, const char *filename ) {

    FILE *fp;
    long i, j, ret;
    char *mxerr0;

    mxerr0 = mxerr1;
    mxerr1 = "MxWrite";

    if(strcmp("stdout", filename)) {
        if( NULL == (fp = fopen(filename, "w")) ) {
            sprintf( mxerr2, "failed open to write \"%s\"", filename);
            PrintError;
            if( NULL != mxerrfoo ) mxerrfoo(); else return EOF;
        }
    } else {
        fp = stdout;
    }

    fprintf(fp, "#  %ld  X  %ld\n", vc, vlen);

    for( i = 0 ; i < vlen ; i++ ) {
        for( j = 0 ; j < vc ; j++ )
            fprintf( fp, mxfmts, *(*(x+j)+i) );
        fprintf( fp, "\n" );
    }

    if(strcmp("stdout", filename)) {
        if( EOF == ( ret = fclose(fp) ) ) {
            sprintf( mxerr2, "failed close upon write on \"%s\"", filename);
            PrintError;
            if( strlen( mxerrs ) ) exit( EXIT_FAILURE );
        }
    } else {
        ret = 0;
    }

    mxerr1 = mxerr0;

    return ret;
}


double **MxAlloc( long vc, long vlen ) {

    long j;
    double **x;
    char *mxerr0;

    mxerr0 = mxerr1;
    mxerr1 = "MxAlloc";

    if( NULL == (x = (double **)calloc(vc, sizeof(double *))) ) {
        sprintf( mxerr2, "no RAM for %ld X %ld", vc, vlen );
        PrintError; TreatError;
    }

    for( j = 0 ; j < vc ; j++ )
        if( NULL == (*(x+j) = (double *)calloc(vlen, sizeof(double))) ) {
            sprintf( mxerr2, "no RAM for %ld X %ld", vc, vlen );
            MxFree(x, vc, j);
            PrintError; TreatError;
        }

    mxerr1 = mxerr0;

    return x;
}


double **MxTrnsp( double **x, long vc, long vlen ) {

    long i, j;
    double **xnew;

    if( NULL != (xnew = MxAlloc(vlen, vc)) )
        for( j = 0 ; j < vc ; j++ )
            for( i = 0 ; i < vlen ; i++ )
                *(*(xnew+i)+j) = *(*(x+j)+i);

    return xnew;
}


double **MxFlttn( double **x, long vc, long vlen ) {

    long i, j;
    double **xnew;

    if( NULL != (xnew = MxAlloc(1, vc*vlen)) )
        for( j = 0 ; j < vc ; j++ )
            for( i = 0 ; i < vlen ; i++ )
                *(*xnew+j*vlen+i) = *(*(x+j)+i);

    return xnew;
}


double **MxMply( double **x, double **y, long vc1, long vlen1, long vlen2 ) {

    long i, j, k;
    double **p;

    if( NULL != (p = MxAlloc(vc1, vlen2)) )
        for( i = 0 ; i < vc1 ; i++ )                            /* rows of x */
            for( k = 0 ; k < vlen2 ; k++ )                   /* columns of y */
                for( j = 0 ; j < vlen1 ; j++ )      /* cols of x = rows of y */
                    *(*(p+i)+k) += *(*(x+i)+j) * *(*(y+j)+k);

    return p;
}


void MxFree( double **x, long vc, long vlen ) {

    long j;

    for( j = 0 ; j < vc ; j++ ) free(*(x+j));  free(x);

    return;
}


void MxCopy( double **from, double **to, long vc, long vlen ) {

    long i, j;

    for( j = 0 ; j < vc ; j++ )
        for( i = 0 ; i < vlen ; i++ )
            *(*(to+j)+i) = *(*(from+j)+i);

}

