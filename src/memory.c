#include <stdlib.h>
#include "../headers/structs.h"

/* Particle Allocator */
struct particle *allocateParticles(int number) {
  struct particle *p = (struct particle *)malloc(number * sizeof(struct particle) );
  return p;
}

/* Grid Allocator */
struct grid *allocateGrid(int numberGridPoints) {
  int i;
  struct grid *g = (struct grid *)malloc( sizeof(struct grid) );

  g->u = (double *)malloc(numberGridPoints * sizeof(double) );
  g->n_i = (double *)malloc(numberGridPoints * sizeof(double) );
  g->n_e = (double *)malloc(numberGridPoints * sizeof(double) );
  g->rho = (double *)malloc(numberGridPoints * sizeof(double) );
  g->J_i = (struct vector2D *)malloc(numberGridPoints * sizeof(struct vector2D) );
  g->J_e = (struct vector2D *)malloc(numberGridPoints * sizeof(struct vector2D) );
  g->J = (struct vector2D *)malloc(numberGridPoints * sizeof(struct vector2D) );

  /* Initialize values */
  for (i=0;i<numberGridPoints;i++) {
    g->u[i] = 0.0;
    g->n_i[i] = 0.0;
    g->n_e[i] = 0.0;
    g->rho[i] = 0.0;
    g->J_i[i].x = 0.0;
    g->J_i[i].y = 0.0;
    g->J_e[i].x = 0.0;
    g->J_e[i].y = 0.0;
    g->J[i].x = 0.0;
    g->J[i].y = 0.0;
  }

  return g;
}

/* Grid De-Allocator */
void deAllocateGrid(struct grid *g) {
  free(g->u); free(g->n_i); free(g->n_e); free(g->rho);
  free(g->J_i); free(g->J_e); free(g->J);
}

/* Field Allocator */
struct field *allocateField(int numberGridPoints) {
  int i; 
  struct field *f = (struct field *)malloc( sizeof(struct field) );

  f->E = (struct vector2D *)malloc(numberGridPoints * sizeof(struct vector2D));
  f->Bz = (double *)malloc(numberGridPoints * sizeof(double));

  /* Initialize values */
  for (i=0;i<numberGridPoints;i++) {
    f->E[i].x = 0.0;
    f->E[i].y = 0.0;
    f->Bz[i] = 0.0;
  }

  return f;
}

/* Field De-Allocator */
void deAllocateField(struct field *f) {
  free(f->E); free(f->Bz);
}
