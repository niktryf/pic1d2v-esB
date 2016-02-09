/*******************************************************************
 *** PIC 1d2v electromagnetic: Poisson 1D Solver
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* A Single 1D Jacobi Iteration - Periodic Boundaries! 
   Solves (d^2/dx^2)u = - rho 
*/
double * jacobiIteration1D (double *u, double *rho, double h, int size)
{
  int i;
  double *u_old;

  /* Allocate u_old array */
  u_old = (double *)malloc((size)*sizeof(double));

  /* Copy u to u_old */
  for(i=0; i<size; i++) u_old[i] = u[i];

  /* Iterate through all inner elements (NOT boundaries) */
  for (i=1; i<size-1; i++) {
    u[i] = 0.5*(u_old[i-1] + u_old[i+1] + h*h*rho[i]);
  }

  /* Clean Up */
  free(u_old);

  return u;
}

/* A Single 1D Gauss-Seidel Iteration - Periodic Boundaries! 
   Solves (d^2/dx^2)u = - rho 
*/
double * gaussSeidelIteration1D (double *u, double *rho, double h, int size)
{
  int i;

  /* Iterate through all inner elements (NOT boundaries) */
  for (i=1; i<size-1; i++) {
    u[i] = 0.5*(u[i-1] + u[i+1] + h*h*rho[i]);
  }

  return u;
}

/* Calculates and returns the residual for the 1D poisson
   problem, as stated in "jacobiIteration1D" function. 
   Uses "quick" jacobi iterative method
   instead of matrix - vector (A * u) multiplication. 
   ----> residual = sqrt ((A*x - r)^2) (CURRENTLY IGNORES BOUNDARIES!)
 */
double residual (double *u, double *rho, double h, int size)
{
  int i;
  double Ax_r, res; 

  res = 0;
  /* Calculate "A*x - rho" and then add to res. */
  for (i=1; i<size-1; i++) {
    Ax_r = (u[i-1] - 2*u[i] + u[i+1] + h*h*rho[i]);
    res += Ax_r*Ax_r;
  }

  return sqrt(res);
}

/* Jacobi Wrapper Function. Normalizes vaccuum permittivity (ε_0) to 1 
   Warning: Assumes boundary conditions have been set!
            This function DOES NOT operate on boundaries!
*/
double * poisson1D (double *u, double *rho, int size, double h)
{
  int i, j, t, maxIterations, nIterations, iterationsPerCheck;
  double res, tolerance;

  /* Set maximum iterations for jacobi solver. 
     We set max = (desired decimal accuracy)/(h^2).
     Since theoretical analysis gives n ~ (2/PI^2)*(d/h^2) ~ 0.2*16/(size^2)
     for the Jacobi method, we set maximum to about 5 
     times more iterations.

     This is probably overkill, and will be reduced when necessary.
  */
  maxIterations = 16*size*size;
  iterationsPerCheck = 100;
  //printf("Maximum Iterations: %d\n", maxIterations);
  
  /* Set tolerance (value of residual that is considered 
     adequate). We set tolerance = size*10^(-d).
   */
  tolerance = size * pow(10.0, -8);
  //printf("Tolerance = %e\n", tolerance);

  nIterations = 0;
  /* Iterate until residual < tolerance or until maxIterations are reached */
    do {
        // Do "iterationsPerCheck" iterations before checking residual
        for (t=0;t<iterationsPerCheck;t++) {
            u = gaussSeidelIteration1D(u, rho, h, size);
            nIterations += 1;
        }
        // Calculate residual
        res = residual(u, rho, h, size);

        // Test print for residual every "iterationsPerCheck"
        //printf("iteration %d: residual = %f\ttolerance: %f\n", nIterations, res, tolerance);

    } while (res > tolerance && nIterations <= maxIterations);


  /* If max iterations are reached, give warning */
  if (nIterations>maxIterations) {
    printf("\n*****************************************\n");
    printf("Warning: Maximum Iterations (%d) reached!\n", maxIterations);
    printf("Residual: %e\n\n", residual(u, rho, h, size));
    printf("*****************************************\n\n");
  }

  return u;
}


