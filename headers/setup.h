struct particle checkPeriodic (struct particle p, double left_bound, double right_bound);

double * applyBoundaryConditions1D (double *a, int size, double leftBound, double rightBound);
double * setBz (double *Bz, int size);

struct particle * setupElectrons(struct particle * p, struct parameters paramT);
struct particle * setupIons(struct particle * p, struct parameters param);
