/********************************************************************
 *** Sets up problem (initializes particles).
 *** The choice to pass the whole parameter structure into the 
 *** setup functions was made in order to be able to work with any 
 *** input parameter. Since these functions are called only once,
 *** at the start of the program, this seems to be a good choice of 
 *** generality vs too much memory moved.
 ***
 ********************************************************************/
#include <stdio.h>
#include <math.h>
#include "../headers/structs.h"
#include "../headers/definitions.h"
#include "../headers/temperature.h"

/*********************************************************************************
 Particle setup
 *********************************************************************************/

/* Applies initial perturbation to ions */
struct particle * perturbIons(struct particle *p, int number, double k) {
  int i;
  double A;

  /* Enter Ion Perturbation here */
  A = 0.5;
  for (i=0;i<number;i++) {
    //p[i].v.x += A*sin(k*2.0*M_PI*i/(number-1));
  }

  return p;
}

/* Applies initial perturbation to electrons */
struct particle * perturbElectrons(struct particle *p, int number, double k) {
  int i;
  double A;

  /* Enter Electron Perturbation here */
  A = 0.5;
  for (i=0;i<number;i++) {
    p[i].v.x += A*sin(k*2.0*M_PI*i/(number-1));
  }

  return p;
}

/* Function that forces periodic conditions:
   Makes particle that goes out of domain 
   appear at the other end.
   ATTENTION: ONLY WORKS IF DISPLACEMENT IS NO MORE THAN ONE GRID LENGTH!!!
   This follows Birdsall and Langdon ES1 approach (Birdsall,Langdon, section 3-7)
 */
struct particle checkPeriodic (struct particle p, double left_bound, double right_bound) 
{
  if (p.r.x > right_bound) {
      p.r.x = p.r.x - (right_bound - left_bound);
    }
    else if (p.r.x < left_bound) {
      p.r.x = p.r.x + (right_bound - left_bound);
    }

  return p;
}

/* Initializes and returns electrons */
struct particle * setupElectrons(struct particle * p, struct parameters param) 
{
  int i;

  /* Apply initial position to electrons */
  for(i=0;i<param.nElectrons;i++) {
    // This sets up electrons uniformly (and NOT on grid points!)
    p[i].r.x = (i+1)*(param.gridEnd - param.gridStart)/(param.nElectrons+1);
    p[i].r.y = 0.0;
  }

  /* Apply initial velocity */
  for(i=0;i<param.nElectrons;i++) {
    p[i].v.x = 0.0;
    p[i].v.y = 0.0;
  }

  /* Apply Temperature */
  p = Maxwell_Boltzmann(p, param.T_e, param.nElectrons, ELECTRON_MASS);

  /* Apply perturbations */
  p = perturbElectrons(p, param.nElectrons, param.k);

  /* Check Periodic Conditions */
  for (i=0; i<param.nElectrons; i++) {
    p[i] = checkPeriodic(p[i], param.gridStart, param.gridEnd);
  }

  return p;
}

/* Initializes and returns ions */
struct particle * setupIons(struct particle * p, struct parameters param) 
{
  int i;

  /* Apply initial position */
  for(i=0;i<param.nIons;i++) {
    // This sets up ions uniformly (and NOT on grid points!)
    p[i].r.x = (i+1)*(param.gridEnd - param.gridStart)/(param.nIons+1);
    p[i].r.y = 0.0;
  }

  /* Apply initial velocity */
  for(i=0;i<param.nIons;i++) {
    p[i].v.x = 0.0;
    p[i].v.y = 0.0;
  }

  /* Apply Temperature (i.e. thermal velocity)*/
  p = Maxwell_Boltzmann(p, param.T_i, param.nIons, ION_MASS);

  /* Apply perturbations */
  p = perturbIons(p, param.nIons, 0.0);

  /* Check Periodic Conditions */
  for (i=0; i<param.nIons; i++) {
    p[i] = checkPeriodic(p[i], param.gridStart, param.gridEnd);
  }

  return p;
}

/*********************************************************************************
 Field setup
 *********************************************************************************/

/* Function that applies 1D boundary conditions to given grid quantity */
double * applyBoundaryConditions1D (double *a, int size, double leftBound, double rightBound) 
{
  a[0] = leftBound; a[size-1] = rightBound;
  return a;
}

/* Sets initial values for Bz (used for uniform magnetic field) */
double * setBz (double *Bz, int size) {

  int i;

  for (i=0;i<size;i++) {
    Bz[i] = 0.0;
  }

  return Bz;
}

