/*******************************************************************
 *** PIC 1d2v electromagnetic: Field functions
 *******************************************************************/
#include <stdlib.h>

#include "../headers/structs.h"

/* Calculates the x component of the Electric Field 
     from the potential (E_x = -du/dx).
   Simple central difference, O(dx^2) error. 
   This is valid only for the x component, for the current model!
   For periodic BCs use periodic treatment of boundaries!
*/
struct vector2D * findEx_fromPotential (struct vector2D *E, double *u, int size, double dx)
{
  int i;

  /* Left Boundary - Downwind scheme (2nd order accuracy) */
  //E[0].x = - (-3*u[0] + 4*u[1] - u[2])/(2*dx);

  /* Left Boundary - Centered Scheme (Periodic!) */
  E[0].x = - (u[1] - u[size-2])/(2*dx);

  /*** Inner points - Central scheme ***/
  for (i=1; i<size-1; i++) {
    E[i].x = - (u[i+1] - u[i-1])/(2*dx);
  }

  /* Right Boundary - Upwind scheme (2nd order accuracy) */
  //E[size-1].x = - (3*u[size-1] - 4*u[size-2] + u[size-3])/(2*dx);

  /* Left Boundary - Centered Scheme (Periodic!) */
  E[size-1].x = - (u[1] - u[size-2])/(2*dx);

  return E;
}
