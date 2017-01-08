//
//  MobilityMatrix.h
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef MobilityMatrix_h
#define MobilityMatrix_h

/* Calculate mobility matrix using RPY tensor */
void regularRPY();

/* Calculate mobility matrix using Ewald summated RPY */
void EwaldSumRPY();

/* Calculate averaged matrix
 * Will also print/update the binary file of the averaged matrix 
 * every 10000 time stpes */
void AvgMatrix(int t);
/* The binary matrix contains a 32 bit (4 bytes) header of number of beads
 * The natrix body contains 3Nc*3Nc single precision floating points (4 bytes each)
 * The whole matrix is read row-wise */

#endif /* MobilityMatrix_h */
