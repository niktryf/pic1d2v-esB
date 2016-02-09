/***
 ***  Structure definitions
 ***/

/* 2D vector structure */
struct vector2D {
  double x;
  double y;
};

/* Particle structure */
struct particle {
  struct vector2D r;
  struct vector2D v;
};

/* grid structure: Holds all quantities that are interpolated
                from the particles to the grid.
*/
struct grid {
  double *u;
  double *n_i, *n_e, *rho;
  struct vector2D *J_i, *J_e, *J;
};

/* field structure: Holds electromagnetic field.
   EM Field is found from grid quantities
   Current program (1d2v): 
   E is 2D (x,y), B is 1D (z)
 */
struct field {
  struct vector2D *E;
  double *Bz;
};

/* parameters structure: Holds everything read from input file 
   (Time, particle and space-grid parameters)

*/
struct parameters {
  double time, dt; 
  int interval;

  int nIons, nElectrons;

  int nGridPoints;
  double gridStart, gridEnd;

  double T_i, T_e, k;
};


