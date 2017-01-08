//
//  MobilityMatrix.c
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "MobilityMatrix.h"
#include "Parameters.h"

void regularRPY(){
    int nn = 3*N;
    
    for(int i = 0; i<nn*nn*NP*NP; ++i){
        matrix->mm[i] = 0.0;
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<NP; ++j){
            int start = i*NP*nn*nn + j*nn*nn;
            for(int ii = 0; ii<N; ++ii){
                for(int jj = 0; jj<N; ++jj){
                    if(i == j && ii == jj){
                        matrix->mm[start+ii*3*nn+jj*3] = 1.0;
                        matrix->mm[start+(ii*3+1)*nn+jj*3+1] = 1.0;
                        matrix->mm[start+(ii*3+2)*nn+jj*3+2] = 1.0;
                    } else{
                        double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                        double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                        double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                        double rr = dx*dx+dy*dy+dz*dz;
                        double r = sqrt(rr);
                        
                        //double mm1,mm2;
                        if(r>=2.0){
                            double mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                            double mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                            matrix->mm[start+ii*3*nn+jj*3] = mm1+mm2*dx*dx;
                            matrix->mm[start+ii*3*nn+jj*3+1] = mm2*dx*dy;
                            matrix->mm[start+ii*3*nn+jj*3+2] = mm2*dx*dz;
                            matrix->mm[start+(ii*3+1)*nn+jj*3] = mm2*dy*dx;
                            matrix->mm[start+(ii*3+1)*nn+jj*3+1] = mm1+mm2*dy*dy;
                            matrix->mm[start+(ii*3+1)*nn+jj*3+2] = mm2*dy*dz;
                            matrix->mm[start+(ii*3+2)*nn+jj*3] = mm2*dz*dx;
                            matrix->mm[start+(ii*3+2)*nn+jj*3+1] = mm2*dz*dy;
                            matrix->mm[start+(ii*3+2)*nn+jj*3+2] = mm1+mm2*dz*dz;
                        }
                        else{
                            double mm1 = 1.0-9.0*r/32.0;
                            double mm2 = 3.0/(32.0*r);
                            matrix->mm[start+ii*3*nn+jj*3] = mm1+mm2*dx*dx;
                            matrix->mm[start+ii*3*nn+jj*3+1] = mm2*dx*dy;
                            matrix->mm[start+ii*3*nn+jj*3+2] = mm2*dx*dz;
                            matrix->mm[start+(ii*3+1)*nn+jj*3] = mm2*dy*dx;
                            matrix->mm[start+(ii*3+1)*nn+jj*3+1] = mm1+mm2*dy*dy;
                            matrix->mm[start+(ii*3+1)*nn+jj*3+2] = mm2*dy*dz;
                            matrix->mm[start+(ii*3+2)*nn+jj*3] = mm2*dz*dx;
                            matrix->mm[start+(ii*3+2)*nn+jj*3+1] = mm2*dz*dy;
                            matrix->mm[start+(ii*3+2)*nn+jj*3+2] = mm1+mm2*dz*dz;
                        }
                        
                        
                    }
                }
            }
        }
    }
}

void AvgMatrix(int t){
    int nn = 3*N;
    
    for(int i = 0; i<nn*nn*NP*NP; ++i){
        matrix->mmAvg[i] += matrix->mm[i];
    }
    
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/MatrixAvg_N%d_hsr%s.bin",directory,N,flowratec);
    
    if(t>0 && t%10000==0){
        FILE *matrixAvg = fopen(str,"wb");
        
        fwrite(&N,sizeof(uint32_t),1,matrixAvg);
    
        for(int i = 0; i<nn*nn*NP*NP; ++i){
            float value = matrix->mmAvg[i]/(t+1);
            fwrite(&value,sizeof(float),1,matrixAvg);
        }
        
        fclose(matrixAvg);
    }
    free(str);
}

void EwaldSumRPY(){
    int nn = 3*N;
    double vol = L*L*L;
    double alpha = 6.0/L;
    
    if(matrix->mm == NULL){
        matrix->mm = calloc(nn*nn*NP*NP, sizeof(float));
    } else{
        int i;
        for(i = 0; i<nn*nn*NP*NP; ++i){
            matrix->mm[i] = 0.0;
        }
    }
    
//#pragma omp parallel for schedule(static)
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<NP; ++j){
            int start = i*NP*nn*nn + j*nn*nn;
            for(int ii = 0; ii<N; ++ii){
                for(int jj = 0; jj<N; ++jj){
                    double dx = Chain[i*N+ii].rx - Chain[j*N+jj].rx;
                    double dy = Chain[i*N+ii].ry - Chain[j*N+jj].ry;
                    double dz = Chain[i*N+ii].rz - Chain[j*N+jj].rz;
                    
                    //Reciprocal Part
                    int kx;
                    for(kx = -KMAX; kx<KMAX+1; ++kx){
                        double rkx = 2.0*M_PI*kx/L;
                        int ky;
                        for(ky = -KMAX; ky<KMAX+1; ++ky){
                            double rky = 2.0*M_PI*ky/L;
                            int kz;
                            for(kz = -KMAX; kz<KMAX+1; ++kz){
                                double rkz = 2.0*M_PI*kz/L;
                                double kk = kx*kx+ky*ky+kz*kz;
                                if(kk != 0){
                                    double mm4 = cos(rkx*dx+rky*dy+rkz*dz);
                                    matrix->mm[start+ii*3*nn+jj*3] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][0]/vol;
                                    matrix->mm[start+ii*3*nn+jj*3+1] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][1]/vol;
                                    matrix->mm[start+ii*3*nn+jj*3+2] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][2]/vol;
                                    matrix->mm[start+(ii*3+1)*nn+jj*3] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][3]/vol;
                                    matrix->mm[start+(ii*3+1)*nn+jj*3+1] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][4]/vol;
                                    matrix->mm[start+(ii*3+1)*nn+jj*3+2] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][5]/vol;
                                    matrix->mm[start+(ii*3+2)*nn+jj*3] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][6]/vol;
                                    matrix->mm[start+(ii*3+2)*nn+jj*3+1] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][7]/vol;
                                    matrix->mm[start+(ii*3+2)*nn+jj*3+2] += mm4*m2[kx+KMAX][ky+KMAX][kz+KMAX][8]/vol;
                                }
                            }
                        }
                    }
                    
                    //Real Part (Cutoff 0.5L)
                    dx -= round(dx/L)*L; dy -= round(dy/L)*L; dz -= round(dz/L)*L;
                    double rr = dx*dx+dy*dy+dz*dz;
                    double r = sqrt(rr);
                    
                    if(i==j && ii == jj){
                        double diag = 1.0-6.0/sqrt(M_PI)*alpha+40.0/3.0/sqrt(M_PI)*alpha*alpha*alpha;
                        matrix->mm[start+ii*3*nn+jj*3] += diag;
                        matrix->mm[start+(ii*3+1)*nn+jj*3+1] += diag;
                        matrix->mm[start+(ii*3+2)*nn+jj*3+2] += diag;
                    }
                    else{
                        double c1 = 0.75/r+0.5*pow(r,-3.0);
                        double c2 = 4.0*pow(alpha,7.0)*pow(r,4.0)+3.0*pow(alpha,3.0)*rr-20.0*pow(alpha,5.0)*rr-4.5*alpha+14.0*pow(alpha,3.0)+alpha/rr;
                        double c3 = 0.75/r-1.5*pow(r,-3.0);
                        double c4 = -4.0*pow(alpha,7.0)*pow(r,4.0)-3.0*pow(alpha,3.0)*rr+16.0*pow(alpha,5.0)*rr+1.5*alpha-2.0*pow(alpha,3.0)-3.0*alpha/rr;
                        
                        double mm1,mm2;
                        if(r>=2.0){
                            mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI);
                            mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI);
                        }
                        else{
                            mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI)+(1.0-9.0*r/32.0-3.0/4.0/r-0.5*pow(r,-3.0));
                            mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI)+(3.0*r/32.0-3.0/4.0/r+1.5*pow(r,-3.0));
                        }
                        matrix->mm[start+ii*3*nn+jj*3] += mm1+mm2*dx*dx/rr;
                        matrix->mm[start+ii*3*nn+jj*3+1] += mm2*dx*dy/rr;
                        matrix->mm[start+ii*3*nn+jj*3+2] += mm2*dx*dz/rr;
                        matrix->mm[start+(ii*3+1)*nn+jj*3] += mm2*dy*dx/rr;
                        matrix->mm[start+(ii*3+1)*nn+jj*3+1] += mm1+mm2*dy*dy/rr;
                        matrix->mm[start+(ii*3+1)*nn+jj*3+2] += mm2*dy*dz/rr;
                        matrix->mm[start+(ii*3+2)*nn+jj*3] += mm2*dz*dx/rr;
                        matrix->mm[start+(ii*3+2)*nn+jj*3+1] += mm2*dz*dy/rr;
                        matrix->mm[start+(ii*3+2)*nn+jj*3+2] += mm1+mm2*dz*dz/rr;
                    }
                }
            }
        }
    }
}