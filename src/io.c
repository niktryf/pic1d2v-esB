/*******************************************************************
 **** PIC - 1d2v electromagnetic - IO Functions
 *******************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include "../headers/structs.h"

#define BUF_LENGTH 100

/* Prints problem parameters to the screen, at the 
   beginning of the program.
 */
void printParameters (struct parameters param, int totalTimeSteps, double dx) 
{
  printf("\n################### PIC Simulation (1d2v) ###################\n");
  printf("# Parameters: \n");
  printf("# \t\tTotal time: %.2f\tdt: %f\n", param.time, param.dt);
  printf("# \t\t(%d timesteps, output every %d steps)\n#\n", totalTimeSteps, param.interval);
  printf("# \t\tNumber of ions: \t%d\n#\t\tNumber of electrons: \t%d\n#\n", param.nIons, param.nElectrons);
  printf("# \t\tIon T: \t\t\t%.3f\n#\t\tElectron T: \t\t%.3f\n#\t\tk: \t\t\t%.2f\n", param.T_i, param.T_e, param.k);
  printf("# \t\tGrid Points: \t\t%d\n#\t\tCell size (dx): \t%f\n#", param.nGridPoints, dx);
  printf("\n#############################################################\n");
}


/* Gets parameters from single line of file input */
struct parameters getParametersFromFile(char * filename) {

  struct parameters p;

  char buf[BUF_LENGTH];

  FILE * inputFile;
  inputFile = fopen(filename, "r");

  while( fgets(buf, BUF_LENGTH, inputFile) != NULL ) {
    /* If scanning Time Parameters */
    if (buf[0] == 'T'){
      sscanf(buf, "%c %lf %lf %d", &buf[0], &p.time, &p.dt, &p.interval);
    }
    /* If scanning Particle Parameters */
    else if (buf[0] == 'P'){
      sscanf(buf, "%c %d %d", &buf[0], &p.nIons, &p.nElectrons);
    }
    /* If scanning Space Parameters */
    else if (buf[0] == 'S'){
      sscanf(buf, "%c %d %lf %lf", &buf[0], &p.nGridPoints, &p.gridStart, &p.gridEnd);
    }
    /* If scanning Other Parameters */
    else if (buf[0] == 'O'){
      sscanf(buf, "%c %lf %lf %lf", &buf[0], &p.T_i, &p.T_e, &p.k);
    }
  }
  
  fclose(inputFile);

  return p;
}

/* 
    Writes scalar quantity for grid points (1D grid)
*/
void writeScalar(double *a, double t, int nGridPoints, char * filename)
{
  int i;
  FILE *outputFile;

  /* Open Output File */
  if (t == 0) {
    outputFile = fopen(filename, "w");
    fprintf(outputFile, "# Gridpoint\tvalue\n\n");
  }
  else outputFile = fopen(filename, "a");
  
  /* Write Output */
  for (i=0; i<nGridPoints; i++) 
  {
    fprintf(outputFile, "%d\t\t%f\n", i, a[i]);
  }
  fprintf(outputFile, "\n");

  /* Close Output File */
  fclose(outputFile);
}

/* 
    Writes 2D Vector Field for grid points (1D grid)
*/
void writeVector(struct vector2D *a, double t, int nGridPoints, char * filename)
{
  int i;
  FILE *outputFile;

  /* Open Output File */
  if (t == 0) {
    outputFile = fopen(filename, "w");
  }
  else outputFile = fopen(filename, "a");
  
  /* Write Output */
  for (i=0; i<nGridPoints; i++) 
  {
    fprintf(outputFile, "%d\t\t%f\t\t%f\n", i, a[i].x, a[i].y);
  }
  fprintf(outputFile, "\n");

  /* Close Output File */
  fclose(outputFile);
}

/* Field output function: prints arrays related to the field (E, B) */
void writeFieldOutput(struct field *f, int nGridPoints, double t) {
  writeScalar(f->Bz, t, nGridPoints, "output/Bz1D.txt");
  writeVector(f->E, t, nGridPoints, "output/E1D.txt");
}

/* Grid output wrapper function: prints arrays related to the grid (rho, u etc) */
void writeGridOutput(struct grid *g, int nGridPoints, double t) {
  writeScalar(g->rho, t, nGridPoints, "output/rho1D.txt");
  writeScalar(g->u, t, nGridPoints, "output/potential1D.txt");
  writeVector(g->J, t, nGridPoints, "output/J1D.txt");
}


