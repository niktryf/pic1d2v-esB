void printParameters (struct parameters param, int totalTimeSteps, double dx);

struct parameters getParametersFromFile(char *filename);

void writeGridOutput(struct grid *g, int nGridPoints, double t);
void writeFieldOutput(struct field *f, int nGridPoints, double t);
