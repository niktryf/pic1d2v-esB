/**************************************************************************
 **** PIC - 1d2v electromagnetic
 **** This file contains functions that interpolate particle data to cells
 **** and grid data to particles.
 **************************************************************************/

#include <stdio.h> //for testing purposes!
#include <math.h>

#include "../headers/structs.h"
#include "../headers/definitions.h"

/****************************************************************
  From grid to particles:
 ****************************************************************/
/* Function that returns the electric field 
   at given particle's position x. It employs CIC
   interpolation.
*/ 
struct vector2D particleE (double x, struct vector2D *E, double dx)
{
  int cell;
  struct vector2D pE;

  /* Find cell index of particle */
  cell = (int)floor(x/dx);

  /* Interpolate fields from neighboring cells */
  pE.x = ( x - cell*dx ) * E[cell+1].x / dx;
  pE.x += ( (cell+1)*dx - x ) * E[cell].x / dx;

  pE.y = ( x - cell*dx ) * E[cell+1].y / dx;
  pE.y += ( (cell+1)*dx - x ) * E[cell].y / dx;

  return pE;
}

/* Function that returns the magnetic field Bz
   at given particle's position x. It employs CIC
   interpolation.
*/ 
double particleBz (double x, double *B, double dx)
{
  int cell;
  double pB;

  /* Find cell index of particle */
  cell = (int)floor(x/dx);

  /* Interpolate fields from neighboring cells */
  pB = ( x - cell*dx ) * B[cell+1] / dx;
  pB += ( (cell+1)*dx - x ) * B[cell] / dx;

  return pB;
}

/* Returns force on given particle.
   Arguments were chosen to satisfy the Runge-Kutta method,
   which relies on incrementing position and velocities on each stage.
*/
struct vector2D particleF(double x, double v_x, double v_y, double particleCharge,
			  struct vector2D *E, double *Bz, double dx)
{
  struct vector2D pE, F;
  double pBz;
  
  /* Interpolate fields to particle */
  pE = particleE(x, E, dx);
  pBz = particleBz(x, Bz, dx);

  /* Calculate force on particle: F = q * [E + (v x B) ]
     The cross product in this 1D version was implemented by hand
     (see the signs). 
  */
  F.x = particleCharge*(pE.x + v_y*pBz);
  F.y = particleCharge*(pE.y - v_x*pBz);

  return F;
}

/****************************************************************
  From particles to grid:
 ****************************************************************/

/* Cloud-In-Cell interpolation (from particle positions to number density).
   First order interpolation - Returns n_(i or e) array
   See Birdsall - Langdon, Part 1, ch.2-6
*/
double * nCIC (double *n, struct particle *p, 
		 double dx, int particleNumber, int nGrid)
{
  int i, cell;

  /* Zero previous calculation */
  for (i=0; i<nGrid; i++) {
    n[i] = 0;
  }

  /* Go through all particles in given species */
  for (i=0; i<particleNumber; i++) {
    //Find cell index of particle
    cell = (int)floor(p[i].r.x/dx);

    /* FOR TESTING - DEBUGGING PURPOSES: */
    //if(cell < 0) printf("\n*** particle %d, cell index < 0 !!! ***\n", i);
    //if(cell >= nGrid) printf("\n*** particle %d, cell index > nGrid !!! ***\n", i);  

    /* Interpolate charge to neighboring cells 
	This is equivalent to counting (but also interpolating)
	At the end n[cell] can be considered as "number of particles in cell"
    */
    n[cell] += ( (cell+1)*dx - p[i].r.x )/dx; 
    n[cell+1] += (p[i].r.x - cell*dx)/dx; 
  }

  /* Calculate actual density by dividing number of particles in cell by dx */
  for(i=0; i<nGrid; i++) {
    n[i] = n[i]/dx;
  }

  /* 
     Force periodic boundaries for particles
     This makes each boundary grid point also take
     into account the other boundary's half
     (Langdon's approach, p. 60/469) 
  */
  n[0] += n[nGrid - 1];
  n[nGrid - 1] = n[0];

  return n;
}

/*
   Finds rho (charge density) from particle positions (interpolation). 
   Also updates n_i, n_e.
   See Birdsall-Langdon Part 1 - Ch.2-6 (p.39/469).
*/
double * interpolateRho (struct grid * g, 
		  struct particle * ions, struct particle * electrons, 
		  struct parameters param, double dx) 
{
  int i;

  /* Interpolation from ions to n_i density */
  g->n_i = nCIC(g->n_i, ions, dx, param.nIons, param.nGridPoints);

  /* Interpolation from electrons to n_e density */
  g->n_e = nCIC(g->n_e, electrons, dx, param.nElectrons, param.nGridPoints);

  /* Calculate charge density */
  for (i=0; i<param.nGridPoints; i++) {
    g->rho[i] = (g->n_i[i]*ION_CHARGE + g->n_e[i]*ELECTRON_CHARGE);
  }

  return g->rho;
}

/******************************************************
 J (current density) interpolation
 ******************************************************/

/* Interpolates current density J for one species */
struct vector2D * jCIC(struct vector2D *j, struct particle *p, 
            double dx, int particleNumber, int nGrid) 
{
  int i, cell;

  /* Zero previous calculation */
  for (i=0; i<nGrid; i++) {
    j[i].x = 0;
    j[i].y = 0;
  }

  /* Go through all particles in given species */
  for (i=0; i<particleNumber; i++) {
    //Find cell index of particle
    cell = (int)floor(p[i].r.x/dx);

    /* TEST - DEBUG - REMOVE */
    //if(cell < 0) printf("\n*** particle %d, charge %f: cell index < 0 !!! ***\n", i, p[i].q);
    //if(cell >= nGrid) printf("\n*** particle %d, charge %f: cell index > nGrid !!! ***\n", i, p[i].q);  

    /* Interpolate charge & velocity to neighboring cells 
	This is equivalent to counting (but also interpolating)
	At the end j[cell] can be considered as "number of particles in cell times velocity"
    */
    j[cell].x += ( (cell+1)*dx - p[i].r.x )*p[i].v.x/dx; 
    j[cell+1].x += (p[i].r.x - cell*dx)*p[i].v.x/dx; 

    j[cell].y += ( (cell+1)*dx - p[i].r.x )*p[i].v.y/dx; 
    j[cell+1].y += (p[i].r.x - cell*dx)*p[i].v.y/dx; 
  }

  /* Force periodic boundaries for particles
     This makes each boundary grid point also take
     into account the other boundary's half
     (Langdon's approach, p. 60/469) */
  j[0].x += j[nGrid - 1].x;
  j[nGrid - 1].x = j[0].x;

  j[0].y += j[nGrid - 1].y;
  j[nGrid - 1].y = j[0].y;

  /* Calculate actual density by dividing number of particles in cell by dx */
  for(i=0; i<nGrid; i++) {
    j[i].x = j[i].x/dx;
    j[i].y = j[i].y/dx;
  }

  return j;
}

/* Calculates Current Density J for both species, by calling jCIC */
struct vector2D * interpolateJ (struct grid * g, struct particle * ions, struct particle * electrons, 
		       struct parameters param, double dx) 
{
  int i;

  /* Zero previous j calculation */
  for (i=0; i<param.nGridPoints; i++) {
    g->J[i].x = 0;
    g->J[i].y = 0;
  }

  /* Interpolation from ions to j_i density */
  g->J_i = jCIC(g->J_i, ions, dx, param.nIons, param.nGridPoints);

  /* Interpolation from electrons to j_e density */
  g->J_e = jCIC(g->J_e, electrons, dx, param.nElectrons, param.nGridPoints);

  /* Calculate charge density */
  for (i=0; i<param.nGridPoints; i++) {
    g->J[i].x = (g->J_i[i].x*ION_CHARGE + g->J_e[i].x*ELECTRON_CHARGE);
    g->J[i].y = (g->J_i[i].y*ION_CHARGE + g->J_e[i].y*ELECTRON_CHARGE);
  }

  return g->J;
}
