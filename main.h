//
//  main.h
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef main_h
#define main_h

#include "Parameters.h"

/* Calculate the nerest image distance of bead i and j */
Vector3D_t getNID(int i, int j);

void applyPBC();
void resetForce();

/* Calculate the boding forces */
void forceSpring();
void forceSpringRing(); //For single ring polymer only

/* Calculate the exclusive volume interactions */
void forceRepulsive();

void updateChainCH();

/* Following functions only used in system with PBC and using TEA for decomposition */
void updateChain();
void locateD_ij(int i, int j, float *D_ij); //pass the selected 3N*3N matrix to *D_ij
void locateForce(int j, float *force); //pass the selected forces to *force

/* Calculate the center of mass of bead i */
Vector3D_t CenterOfMass(int i);

/* Store the center of mass of all the beads in the system in *R_CoM */
void storeCoM();

/* Following functions calculate varies of properties
 * I comment the print functions due to the file naming systmem I was using
 * Just create ur own print functions */
void LongTimeDc(int t); //long-time limit self diffusivity
void Rgyration(int t); //Calculate radii of gyration and end-to-end distances
void Rendtoend(int t); //Calculate end-to-end vectors
void Extension(int t); //Calcualate extension (in x direction only in this fcn)
void ExtensionRing(int t); //Calculate extension of a ring polymer (in all diretions)

float gasdev(long *idum);

void printTrajectory(int t);

/* Track wall-clock time, output in unit of nano second */
long long timer();

#endif /* main_h */
