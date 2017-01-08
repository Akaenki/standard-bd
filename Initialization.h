//
//  Initialization.h
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Initialization_h
#define Initialization_h

void initialization(int argc, char * argv[]);

/* Function used to read intial traj from a file */
void readChains();

/* Generate a .txt output file for record */
void generateOutput();

long initRan();

/* Initialize the parallel environment */
void initParaEvir(int numThreads);

void initChain();

/* Initial ring structure, add it after the function initChain()
 * Only works for single ring in this file */
void initRingStructure();

/* Initialize the m2 matrix in Ewald sum calculations */
void initM2();

float ran1(long *idum);

void printInitial();

#endif /* Initialization_h */
