/*******************************************************************
 *** PIC 1d2v electromagnetic: Particle Mover (Runge Kutta 4)
 *******************************************************************/

#include "../headers/structs.h"
#include "../headers/interpolate.h"

/* Returns the right-hand side of the particle equation of motion:
   rhs -> (Force vector)/mass.
 */
struct vector2D particleRHS(double x, double v_x, double v_y, 
			    double particleCharge, double particleMass,
			    struct vector2D *E, double *Bz,
			    double dx)
{
  struct vector2D F, rhs;

  F = particleF(x, v_x, v_y, particleCharge, E, Bz, dx);

  rhs.x = F.x/particleMass;
  rhs.y = F.y/particleMass;

  return rhs;
}

/* 
   Basic Runge Kutta 4 function: 
   Returns new particle from old 
   (position and velocity in vector2D form).
   Specifically made for PIC 1d2v model ( F = F(x, v_x, v_y, E, B) ).
*/
struct particle moveParticle(struct particle p, double charge, double mass, 
			     struct field *f, 
			     double dx, double h)
{
  struct vector2D k[4], l[4], rhs;
  
  // Stage 1
  rhs = particleRHS(p.r.x, p.v.x, p.v.y, charge, mass, f->E, f->Bz, dx);
  k[0].x = h*rhs.x; //velocity
  k[0].y = h*rhs.y; //velocity
  l[0].x = h*p.v.x; //position
  l[0].y = h*p.v.y; //position

  // Stage 2
  rhs = particleRHS(p.r.x + 0.5*l[0].x, p.v.x + 0.5*k[0].x, p.v.y + 0.5*k[0].y, 
		    charge, mass, f->E, f->Bz, dx);
  k[1].x = h*rhs.x; 
  k[1].y = h*rhs.y;
  l[1].x = h*(p.v.x + 0.5*k[0].x); 
  l[1].y = h*(p.v.y + 0.5*k[0].y);

  // Stage 3
  rhs = particleRHS(p.r.x + 0.5*l[1].x, p.v.x + 0.5*k[1].x, p.v.y + 0.5*k[1].y, 
		    charge, mass, f->E, f->Bz, dx);
  k[2].x = h*rhs.x; 
  k[2].y = h*rhs.y;
  l[2].x = h*(p.v.x + 0.5*k[1].x); 
  l[2].y = h*(p.v.y + 0.5*k[1].y);

  // Stage 4
  rhs = particleRHS(p.r.x + l[2].x, p.v.x + k[2].x, p.v.y + k[2].y, 
		    charge, mass, f->E, f->Bz, dx);
  k[3].x = h*rhs.x; 
  k[3].y = h*rhs.y;
  l[3].x = h*(p.v.x + k[2].x); 
  l[3].y = h*(p.v.y + k[2].y);

  // Calculate new r, v:
  p.v.x += (1.0/6.0)*( k[0].x + 2*(k[1].x + k[2].x) + k[3].x );
  p.r.x += (1.0/6.0)*( l[0].x + 2*(l[1].x + l[2].x) + l[3].x );

  p.v.y += (1.0/6.0)*( k[0].y + 2*(k[1].y + k[2].y) + k[3].y );
  p.r.y += (1.0/6.0)*( l[0].y + 2*(l[1].y + l[2].y) + l[3].y );

  return p;
}










