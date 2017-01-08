//
//  Parameters.h
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* This file contains all the folloqings
 * 1. common headers
 * 2. stucture definations
 * 3. all the global variables */

#ifndef Parameters_h
#define Parameters_h

#if defined(HYBRID_PARALLEL)
#include <mpi.h>
#endif

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cblas.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define EP 0.9624
#define THETA0 31.7

#ifndef M_PI
#define M_PI 3.1415926535
#endif


/////////////////////////////////////////////////////////////////////////////
// Structures
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    double x,y,z;
} Vector3D_t;

typedef struct {
    double rx,ry,rz; //beads location
    double fx,fy,fz; //forces acting on beads
} Chain_t;

/* Mobility Matrix */
typedef struct {
    float* mm; //Matrix updated
    float* mmAvg; //Matrix averaged
} Matrix_t;

/* Geyer-Winter(TEA) Parameters */
typedef struct {
    double beta_ij;
    float* C_i;
} gwParm_t;

/////////////////////////////////////////////////////////////////////////////
// Variables/Arrays
/////////////////////////////////////////////////////////////////////////////
long *idum; //SEED

int flowType; //flow type: 0 - static, 1 - elongational, 2 - steady shear
double flowRate;

long Starttime;

uint32_t N;
int numThread;

Vector3D_t *R_CoM; //stores center of mass of all the beads
double MSD_CoM; //mean square distance of COM
Vector3D_t ext; //extension 

Chain_t* Chain;
Matrix_t* matrix;

gwParm_t* gwParm;

float *lm; //Cholesky Decomposition
float ****m2; //the pre-defined m2 for ewald sum calculation


/////////////////////////////////////////////////////////////////////////////
// Others
/////////////////////////////////////////////////////////////////////////////

int trajStep;

double c_star; //Overlap concentration (normalized)
double c_normal; //normalized concentration

double Rg,Ree;

char* outputName;
char* trajName;
char* directory;
char* flowratec;

bool ReadInitFromFile;
bool CheckOverlap;
bool GeyerWinter;
bool usePBC;

#endif /* Parameters_h */
