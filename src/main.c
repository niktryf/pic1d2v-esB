#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "../headers/structs.h"
#include "../headers/memory.h"
#include "../headers/io.h"
#include "../headers/setup.h"
#include "../headers/wrappers.h"
#include "../headers/fields.h"
#include "../headers/mover.h"

#include "../headers/definitions.h"

int main() {

  int i, t, output;
  double tStart, tEnd;

  /* Get parameters from input file */
  struct parameters param;
  param = getParametersFromFile("input.txt");

  /* Calculate secondary parameters */
  int totalTimeSteps = (int)(param.time/param.dt);
  int nOutput = totalTimeSteps/param.interval;
    if (totalTimeSteps%param.interval!=0) nOutput +=1;
  double dx = (param.gridEnd - param.gridStart)/(param.nGridPoints - 1);

  /* Print Simulation Parameters */
  printParameters (param, totalTimeSteps, dx);

  /*** Memory Allocation ***********************/
  /* Declare and Allocate Memory for Particles */
  struct particle *ions, *electrons;
  ions = allocateParticles(param.nIons);
  electrons = allocateParticles(param.nElectrons); 

  /* Declare and Allocate Memory for Grid */
  struct grid *g;
  g = allocateGrid(param.nGridPoints);

  /* Declare and Allocate Memory for Field */
  struct field *f;
  f = allocateField(param.nGridPoints);
  /*** Memory Allocation End *******************/


  /************* Setup ********************/
  /* Particles: */
  electrons = setupElectrons(electrons, param);
  ions = setupIons(ions, param);

  /* Apply boundary conditions (potential): */
  g->u = applyBoundaryConditions1D (g->u, param.nGridPoints, 0.0, 0.0);

  /* Set initial fields */
  // Set Steady/uniform Bz
  f->Bz = setBz(f->Bz, param.nGridPoints);

  /* Start timing */
  tStart = omp_get_wtime();

  /**** START ITERATING ****/   
  for (output=0;output<=nOutput;output++) { 
    /* Write output */
    writeGridOutput(g, param.nGridPoints, output);
    writeFieldOutput(f, param.nGridPoints, output);   
    for (t=0; t<param.interval; t++) {
      
      /* Calculate grid quantities (interpolate n_i, n_e -> rho, j_i, j_e -> J and solve for potential u) */
      g = fromParticlesToGrid(g, ions, electrons, param, dx);

      /* Differentiate potential (u) to get the Electric Field E_x ( du/dx = -E(x) )*/
      f->E = findEx_fromPotential (f->E, g->u, param.nGridPoints, dx); 

      /* Move ions and electrons with the new values for E_x */
      for(i=0;i<param.nIons;i++) {
        ions[i] = moveParticle(ions[i], ION_CHARGE, ION_MASS, f, dx, param.dt);
        ions[i] = checkPeriodic (ions[i], param.gridStart, param.gridEnd);
      }
      for(i=0;i<param.nElectrons;i++) {
        electrons[i] = moveParticle(electrons[i], ELECTRON_CHARGE, ELECTRON_MASS, f, dx, param.dt);
        electrons[i] = checkPeriodic (electrons[i], param.gridStart, param.gridEnd);
      }
    }
  }

  /* Stop timing */
  tEnd = omp_get_wtime();

  /* End */
  printf("\n...done! Time: %f sec. \n\n", tEnd - tStart);

  /*** Free memory ***/
  free(ions); free(electrons); 
  deAllocateGrid(g); free(g);
  deAllocateField(f); free(f);

  return 0;
}
