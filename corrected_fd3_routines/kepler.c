
#include "kepler.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define OVERPHASE(x)  ( 2.0 * M_PI * floor ( 0.5 + 0.5 * x * M_1_PI ) )

#define DEBUG 0

/****************************************************************************/

struct delta_params { double ecosmu, esinmu; };

static double delta_f ( double delta, void *params ) {

  struct delta_params * op = ( struct delta_params * ) params;
  double ecosmu = op->ecosmu;
  double esinmu = op->esinmu;

  return ecosmu * sin(delta) + esinmu * cos(delta) - delta;

}

static double delta_df ( double delta, void *params ) {

  struct delta_params * op = ( struct delta_params * ) params;
  double ecosmu = op->ecosmu;
  double esinmu = op->esinmu;

  return ecosmu * cos(delta) - esinmu * sin(delta) - 1.0;

}

static void delta_fdf ( double delta, void *params, double * f, double * df ) {
  struct delta_params * op = ( struct delta_params * ) params;
  double ecosmu = op->ecosmu;
  double esinmu = op->esinmu;
  double cosdelta = cos(delta);
  double sindelta = sin(delta);
  *f  = ecosmu * sindelta + esinmu * cosdelta - delta;
  *df = ecosmu * cosdelta - esinmu * sindelta - 1.0;
}

double kepler_psiofmu ( double mu, double ecc ) {

  int status;
  int iter = 0, max_iter = 100;
  double x0, x, x_expected;
  gsl_function_fdf FDF;
  static gsl_root_fdfsolver * solver;
  static int solverexists = 0;
  struct delta_params params = { ecc * cos(mu), ecc * sin(mu) };

  if ( ! solverexists ) { solverexists = 1;
    solver = gsl_root_fdfsolver_alloc ( gsl_root_fdfsolver_newton );
  }

  if ( DEBUG ) printf( "kepler_psiofmu ( mu = %lg, ecc = %lg )\n", mu, ecc );

  /* early exit make this more tollerant */
  if ( 0.0 == ecc ) return mu;
  if ( 0.0 == mu - OVERPHASE(mu) ) return 0.0;

  FDF.f = &delta_f;
  FDF.df = &delta_df;
  FDF.fdf = &delta_fdf;
  FDF.params = &params;

  x_expected = x = (&params)->esinmu;

  gsl_root_fdfsolver_set (solver, &FDF, x);

  if ( DEBUG ) printf ("using %s method\n", gsl_root_fdfsolver_name (solver));
  if ( DEBUG ) printf ("%-5s %10s %10s\n", "iter", "root", "err");

  do {
    iter++;
    status = gsl_root_fdfsolver_iterate (solver);
    x0 = x;
    x = gsl_root_fdfsolver_root (solver);
    status = gsl_root_test_delta (x, x0, 0, 1e-3);
    if ( DEBUG && status == GSL_SUCCESS ) printf ("Converged:\n");
    if ( DEBUG ) printf ("%5d %10.7f %+10.7f\n", iter, x, x - x_expected);
  } while (status == GSL_CONTINUE && iter < max_iter);

  if ( DEBUG ) printf( "kepler_psiofmu ( mu = %lg, ecc = %lg ) = %lg\n",
    mu, ecc, mu+x );

  return  mu + x;

}

/****************************************************************************/

double kepler_thetaofmu ( double mu, double ecc ) {

  double overmu = OVERPHASE(mu), psi, theta;

  psi = kepler_psiofmu ( mu, ecc );
  theta = 2.0 * atan( sqrt((1.0+ecc)/(1.0-ecc)) * tan(0.5*psi) );

  return  overmu + theta;

}

/****************************************************************************/

double kepler_rv ( double mu, double ecc, double omega ) {

  double psi, theta;

  psi = kepler_psiofmu ( mu, ecc );
  theta = 2.0 * atan( sqrt((1.0+ecc)/(1.0-ecc)) * tan(0.5*psi) );

  return  cos(theta + omega) + ecc*cos(omega);
}

/****************************************************************************/

double kepler_lite ( double mu, double ecc, double omega ) {

  double psi, theta;

  psi = kepler_psiofmu ( mu, ecc );
  theta = 2.0 * atan( sqrt((1.0+ecc)/(1.0-ecc)) * tan(0.5*psi) );

  return  (1.0 - ecc*ecc) * sin(theta + omega) / (1.0 + ecc*cos(theta));
}

/****************************************************************************/

double kepler_sep ( double mu, double ecc, double omega, double inc ) {

  double psi, theta, tmpd;

  psi = kepler_psiofmu ( mu, ecc );
  theta = 2.0 * atan( sqrt((1.0+ecc)/(1.0-ecc)) * tan(0.5*psi) );
  tmpd = sin(theta + omega) * sin(inc);
  tmpd = sqrt(1.0 - tmpd*tmpd);

  return  tmpd * (1.0 - ecc*ecc) / (1.0 + ecc * cos(theta));
}

/****************************************************************************/

