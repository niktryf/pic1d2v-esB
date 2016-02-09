/**********************************************************
 ***
 *** Wrapper functions for PIC - 1d2v 
 ***
 **********************************************************/

#include "../headers/structs.h"
#include "../headers/poisson.h"
#include "../headers/interpolate.h"

/* 
  Evaluates grid quantities from particles
*/
struct grid * fromParticlesToGrid(struct grid *g, struct particle *ions, struct particle *electrons, 
                                             struct parameters param, double dx) 
{
  /* Charge (rho) interpolation, from particles to grid */
  g->rho = interpolateRho (g, ions, electrons, param, dx);

  /* Solution of Poisson Equation: rho -> u -> E_x */
  g->u = poisson1D (g->u, g->rho, param.nGridPoints, dx);

  /* Current density (J) interpolation, from particles to grid */
  g->J = interpolateJ (g, ions, electrons, param, dx);

  return g;
}
