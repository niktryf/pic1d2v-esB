/*** Header files for functions in memory.c ***/

struct particle *allocateParticles(int number);

struct grid *allocateGrid(int numberGridPoints);
void deAllocateGrid(struct grid *g);

struct field *allocateField(int numberGridPoints);
void deAllocateField(struct field *f);

