
#include <stdlib.h>
#include <math.h>

#include "kepler.h"
#include "triorb.h"

#define LITESPEED 2.998E05

int triorb_rv ( double *op, double t, double *rv ) {

  double wp, wt, we, woAB, waAB, waC;
  double tp, tt, te, toA, tkA, tkB, tdo;
  double wliteAB, wrvAB, trvA, L, dLdt, ct, rv_C, rv_AB, rv_A, rv_B;

  /* Wide orbit */

  wp   = op[0];
  wt   = op[1];
  we   = op[2];
  woAB = op[3];
  waAB = op[4];
  waC  = op[5];

  wliteAB = kepler_lite ( 2 * M_PI * ( t - wt ) / wp, we, woAB );
  wrvAB = kepler_rv ( 2 * M_PI * ( t - wt ) / wp, we, woAB );

  /* tight orbit */

  tp   = op[6];
  tt   = op[7];
  te   = op[8];
  toA  = op[9];
  tkA  = op[10];
  tkB  = op[11];
  tdo  = op[12];

  /* C */

  L = - waC * wliteAB;
  dLdt = - 2 * M_PI * waC * wrvAB / ( wp * sqrt ( 1 - we*we ) );
  ct = (t - L + dLdt * t) / (1 + dLdt);
  rv_C = LITESPEED * 2 * M_PI
    * waC * kepler_rv ( 2 * M_PI * ( ct - wt ) / wp, we, M_PI + woAB )
     / ( wp * sqrt ( 1 - we*we ) );

  /* AB */

  L = waC * wliteAB;
  dLdt = 2 * M_PI * waAB * wrvAB / ( wp * sqrt ( 1 - we*we ) );
  ct = (t - L + dLdt * t) / (1 + dLdt);
  rv_AB = LITESPEED * 2 * M_PI
    * waAB * kepler_rv ( 2 * M_PI * ( ct - wt ) / wp, we, woAB )
      / ( wp * sqrt ( 1 - we*we ) );

  /* tight orbit */

  trvA = kepler_rv ( 2 * M_PI * ( ct - tt ) / tp, te, toA );

  rv_A = rv_AB + tkA * trvA;
  rv_B = rv_AB - tkB * trvA;

  /* finish */

  *(rv+0) = rv_A;
  *(rv+1) = rv_B;
  *(rv+2) = rv_C;

  return EXIT_SUCCESS;
}

