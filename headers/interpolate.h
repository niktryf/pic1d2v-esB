double * interpolateRho (struct grid * g, 
		  struct particle * ions, struct particle * electrons, 
		  struct parameters param, double dx);

struct vector2D * interpolateJ (struct grid * g, struct particle * ions, struct particle * electrons, 
		       struct parameters param, double dx);

struct vector2D particleF(double x, double v_x, double v_y, double particleCharge,
			  struct vector2D *E, double *Bz, double dx);
