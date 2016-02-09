/**********************************************************
 ***
 *** Maxwell-Boltzmann thermal velocity for PIC - 1d2v 
 ***
 **********************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../headers/structs.h"

#define K_B 1.0

/*********************************************************************
 Maxwell - Boltzmann Velocity Distribution Functions
 *********************************************************************/
/* Returns Random Uniform between 0 and 1 */
double uniformRandom()
{
  return ( ((double) rand()) / RAND_MAX );
}

/* Box Muller Method: Returns Random Number from Normal Distribution */
double BoxMuller()
{
  double u1=uniformRandom();
  double u2=uniformRandom();
  return sqrt(-2.0*log(u1))*cos(2*M_PI*u2); 
}

/* Apply Desired Average and Standard Deviation 
   To Given Array (of structures) With Values From Normal Distribution
   (-> Initial Average = 0, Initial Standard Deviation = 1)
 */
struct particle * Average_StdDev(struct particle *p, double average, double stddev, int number) 
{
  int i;

  /* Apply Standard Deviation First */
  for (i=0; i<number; i++) {
    p[i].v.x *= sqrt(stddev);
  }

  /* Apply Average Second */
  for (i=0; i<number; i++) {
    p[i].v.x += average;
  }

  return p;
}

/* Calculates average from array */
double average (struct particle *p, int size) 
{
  int i;
  double av = 0;

  for (i=0; i<size; i++) {
    av += p[i].v.x;
  }

  av/=size;

  return av;
}

/* Calculates Standard Deviation from array, given average */
double standardDev(struct particle *p, double average, int size)
{
  int i;
  double stddev = 0;

  for (i=0; i<size; i++) {
    stddev += (average - p[i].v.x)*(average - p[i].v.x);
  }

  stddev = stddev / size;

  return stddev;
}

/* Sets Maxwell-Boltzmann velocities to array of particles of given mass (struct particle) */
struct particle * Maxwell_Boltzmann(struct particle *p, double T, int number, double mass) {
  
  int i;

  /* Seed RNG */
  //srand(time(NULL));
  srand(101);
  
  /* Fill array with Gaussian Distribution Values 
    (average: 0, standard deviation: 1)
  */
  for (i=0; i<number; i++) {
    p[i].v.x = BoxMuller();
  }

  /* Apply Desired Average and Standard Deviation for M_B distribution:
     average: desired (if any) beam velocity
     std_dev = sqrt(kT/m)
     *** k (Boltzmann's constant) set to 1.0 here! ***
  */
  p = Average_StdDev(p, 0, sqrt( (K_B*T)/mass ), number);

  
  /*** Below: Testing Purposes for M-B distribution - May Remove Later ***/
  /* Calculate average */
  double av = average(p, number);

  /* Calculate Standard Deviation */
  double stddev = standardDev(p, av, number);

  /* Print Maxwell-Boltzman Test Results 
  printf("\nMaxwell - Boltzmann Distribution: \n");
  printf("Average = %f\n", av);
  printf("Standard Deviation = %f\n\n", stddev); */

  return p;
}
/*********************************************************************
 Maxwell - Boltzmann Velocity Distribution Functions (END)
 *********************************************************************/
