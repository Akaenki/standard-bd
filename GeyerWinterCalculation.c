//
//  GeyerWinterCalculation.c
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* this file contains functions calculate the decomposition of mobility matrix 
 * Not only the gw calculation (even it's named so)*/

#include "GeyerWinterCalculation.h"
#include "Parameters.h"

void GeyerWinterCalculation(){
    double eps = 0.0;
    int nn = 3*N;
    int nt = nn*NP; /*dimension*/
    
    int i;
    for(i = 0; i<nt*nt; ++i){
        eps += matrix->mm[i];
    }
    
    eps -= nt; /*exclude diagonal(self) components*/
    
    eps /= nt*nt - nt;
    
    gwParm->beta_ij = (1-sqrt(1-((nt-1)*eps*eps-(nt-2)*eps)))/((nt-1)*eps*eps-(nt-2)*eps);
    //printf("beta=%f\n", beta_ij);
    
    for(int i = 0; i<nt; ++i){
        gwParm->C_i[i] = 0.0;
    }
    
    for(i = 0; i<nt; ++i){
        int j;
        for(j = 0; j<nt; ++j){
            gwParm->C_i[i] += matrix->mm[i*nt+j]*matrix->mm[i*nt+j];
        }
    }
    
    for(i = 0; i<nt; ++i){
        gwParm->C_i[i] -= 1.0;
        gwParm->C_i[i] = gwParm->C_i[i]*gwParm->beta_ij*gwParm->beta_ij+1.0;
        gwParm->C_i[i] = 1.0/sqrt(gwParm->C_i[i]);
        //printf("C%d=%f\n",i,gwParm->C_i[i]);
    }
}

void Cholesky(int t){
    int n = 3*N*NP;
    int error = 0;
    if(lm == NULL){
        lm = calloc(n*n,sizeof(float));
    } else{
        int i;
        for(i = 0; i<n*n; ++i){
            lm[i] = 0.0;
        }
    }
    
    int i;
    for (i = 0; i<n; i++){
        int j;
        for (j = 0; j<(i+1); j++){
            double s = 0;
            int k;
            for (k = 0; k < j; k++){
                s += lm[i*n+k]*lm[j*n+k];
            }
            
            if(i == j){
                if((matrix->mm[i*n+i]-s)<=0){
                    error = 1;
                }
                lm[i*n+j] = sqrt(matrix->mm[i*n+i]-s);
            } else{
                lm[i*n+j] = (matrix->mm[i*n+j]-s)/lm[j*n+j];
            }
        }
    }
    
    if(error){
        printf("t=%d\n",t);
        for(i = 0; i<n; ++i){
            int j;
            for(j = 0; j<n; ++j){
                lm[i*n+j] = 0.0;
                if(i == j) lm[i*n+j] = 1.0;
            }
        }
    }
}