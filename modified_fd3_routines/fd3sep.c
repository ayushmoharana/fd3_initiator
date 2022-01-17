
#include <math.h>
#include <stdlib.h>

#include "fd3sep.h"
/* #include "mxfuns.h" */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#define SVCUT 1.0e-9

/*************************************************************************/

double fd3sep ( long K, long M, long N, double **dftobs, double *sig,
   double **rvm, double **lfm, double **dftmod, double **dftres ) {

   long i, j, k, n;
   double s2;
   gsl_matrix *A, *U, *X, *V;
   gsl_vector *S, *w, *b, *x;

   /* prepare */

   A = gsl_matrix_alloc ( 2*M, 2*K );
   U = gsl_matrix_alloc ( 2*M, 2*K );
   X = gsl_matrix_alloc ( 2*K, 2*K );
   V = gsl_matrix_alloc ( 2*K, 2*K );
   S = gsl_vector_alloc ( 2*K );
   w = gsl_vector_alloc ( 2*K );
   b = gsl_vector_alloc ( 2*M );
   x = gsl_vector_alloc ( 2*K );

   /* for each DFT component */

   for ( s2 = n = 0 ; n <= N/2 ; n++ ) {

      /* assemble data vector and model matrix */

      for ( j = 0 ; j < M ; j++ ) {
         double s = *(sig+j);
         gsl_vector_set ( b, 2*j,   *(*(dftobs+j)+2*n)  /s );
         gsl_vector_set ( b, 2*j+1, *(*(dftobs+j)+2*n+1)/s );
         for ( k = 0 ; k < K ; k++ ) {
            double q, v, fv, rez, imz;
            q = 2.0 * M_PI * ((double)n)/((double)N);
            v = *(*(rvm+k)+j);
            fv = floor(v);
            rez = *(*(lfm+k)+j) *
                ( ((fv+1)-v)*cos(fv*q) + (v-fv)*cos((fv+1)*q) );
            imz = - *(*(lfm+k)+j) *
                ( ((fv+1)-v)*sin(fv*q) + (v-fv)*sin((fv+1)*q) );
            gsl_matrix_set ( A, 2*j,   2*k,     rez/s );
            gsl_matrix_set ( A, 2*j,   2*k+1, - imz/s );
            gsl_matrix_set ( A, 2*j+1, 2*k,     imz/s );
            gsl_matrix_set ( A, 2*j+1, 2*k+1,   rez/s );
         }
         gsl_matrix_memcpy ( U, A );
      }

      /* solve for model parameters */

      /* gsl_linalg_SV_decomp ( U, V, S, w ); */
      /* gsl_linalg_SV_decomp_mod ( U, X, V, S, w ); */
      /* gsl_linalg_SV_decomp_jacobi ( U, V, S ); */

      gsl_linalg_SV_decomp_mod ( U, X, V, S, w );
      for ( k = 0 ; k < 2*K-1 ; k++ )
         if ( gsl_vector_get(S,2*K-1-k)/gsl_vector_get(S,0) < SVCUT )
            gsl_vector_set ( S, 2*K-1-k, 0 );
      gsl_linalg_SV_solve ( U, V, S, b, x );

      /* copy to output */

      for ( k = 0 ; k < K ; k++ ) {
         *(*(dftmod+k)+2*n  ) = gsl_vector_get ( x, 2*k   );
         *(*(dftmod+k)+2*n+1) = gsl_vector_get ( x, 2*k+1 );
      }

      /* compute s2 */

      for ( j = 0 ; j < M ; j++ ) for ( i = 0 ; i < 2 ; i++ ) {

        double bc, db, s = *(sig+j);

        for ( bc = k = 0 ; k < K ; k++ ) {
          bc += gsl_matrix_get( A, 2*j+i, 2*k   ) * gsl_vector_get( x, 2*k   );
          bc += gsl_matrix_get( A, 2*j+i, 2*k+1 ) * gsl_vector_get( x, 2*k+1 );
        }
        db = gsl_vector_get ( b, 2*j+i ) - bc;
        s2 += db * db * ( n % ((N+1)/2) ? 2 : 1 );
        *(*(dftres+j)+2*n+i) = db * s;

      }

   }

   /* close */

   gsl_matrix_free ( A );
   gsl_matrix_free ( U );
   gsl_matrix_free ( X );
   gsl_matrix_free ( V );
   gsl_vector_free ( S );
   gsl_vector_free ( w );
   gsl_vector_free ( b );
   gsl_vector_free ( x );

   return s2;

}

/*************************************************************************/

void dft_fwd ( long m, long n, double **mxin, double **mxout ) {

   long i, j;
   double a = 1.0 / sqrt(n);
   gsl_fft_real_wavetable * dft_rewt;
   gsl_fft_real_workspace * dft_rews;

   dft_rewt = gsl_fft_real_wavetable_alloc (n);
   dft_rews = gsl_fft_real_workspace_alloc (n);

   for ( j = 0 ; j < m ; j++ ) {
      for ( i = 0; i < n; i++ ) *(*(mxout+j)+i+1) = a * *(*(mxin+j)+i);
      gsl_fft_real_transform ( *(mxout+j)+1, 1, n, dft_rewt, dft_rews );
      *(*(mxout+j)+0) = *(*(mxout+j)+1);
      *(*(mxout+j)+1) = 0;
      if ( ! n % 2 ) *(*(mxout+j)+n+1) = 0; /* if n even */
   }

   gsl_fft_real_wavetable_free ( dft_rewt );
   gsl_fft_real_workspace_free ( dft_rews );

}

/*************************************************************************/

void dft_bck ( long m, long n, double **mxin, double **mxout ) {

   long i, j;
   double a = 1.0 / sqrt(n);
   gsl_fft_halfcomplex_wavetable * dft_hcwt;
   gsl_fft_real_workspace * dft_rews;

   dft_hcwt = gsl_fft_halfcomplex_wavetable_alloc (n);
   dft_rews = gsl_fft_real_workspace_alloc (n);

   for ( j = 0 ; j < m ; j++ ) {
      *(*(mxout+j)+0) = *(*(mxin+j)+0) / a;
      for ( i = 1; i < n; i++ ) *(*(mxout+j)+i) = *(*(mxin+j)+i+1) / a;
      gsl_fft_halfcomplex_inverse ( *(mxout+j), 1, n, dft_hcwt, dft_rews );
   }

   gsl_fft_halfcomplex_wavetable_free ( dft_hcwt );
   gsl_fft_real_workspace_free ( dft_rews );

}

/*************************************************************************/

