//
//  Initialization.c
//  EwaldSumRPY
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* this file contains functions that initialize the simulation */

#include "Initialization.h"
#include "Parameters.h"
#include "util.h"

void initialization(int argc, char * argv[]){
    ParseInput(argc, argv);
    
    trajStep = 1000; //frequency of printing the trajectories
    
    GeyerWinter = true; //where you want to use
    
    directory = "./output"; //path of output files
    
    usePBC = false; //change it to true for smdlt systems
    
    //ReadInitFromFile = false; //if you want to read the initial traj from a file
    
    CheckOverlap = true; //can turn off it when system is too dense
    
    matrix = (Matrix_t*)malloc(sizeof(Matrix_t));
    matrix->mm = (float*)malloc(9*N*N*NP*NP*sizeof(float));
    matrix->mmAvg = (float*)calloc(9*N*N*NP*NP,sizeof(float));
    gwParm = (gwParm_t*)malloc(sizeof(gwParm_t));
    gwParm->C_i = (float*)malloc(3*N*NP*sizeof(float));
    
    outputName = malloc(100*sizeof(char));
    sprintf(outputName,"%s/sc_N%d_hsr%s.txt",directory,N,flowratec);
    trajName = malloc(100*sizeof(char));
    sprintf(trajName,"%s/sc_N%d_hsr%s.xyz",directory,N,flowratec);
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    idum = malloc(sizeof(long));
    *idum = -1;
    ran1(idum);
    *idum = -1*initRan();
    
    initParaEvir(numThread);
    
    initChain();
    //if(ReadInitFromFile)
    //readChains();
    
    if(usePBC) initM2();
    
    MSD_CoM = 0.0;
    
    ext.x = 0.0; ext.y = 0.0; ext.z = 0.0;
    
    //generateOutput();
    
    //printInitial();
}

void readChains(){
    Chain = calloc(N*NP, sizeof(Chain_t));
    
    char* str = malloc(100*sizeof(char));
    sprintf(str,"coil_%d.xyz",N);
    
    FILE *trajec;
    trajec = fopen(str,"r");
    for(int i = 0; i<N*NP; ++i){
        fscanf(trajec, "A %lf %lf %lf\n", &Chain[i].rx, &Chain[i].ry, &Chain[i].rz);
    }
    fclose(trajec);
    
    free(str);
}

void generateOutput(){
    FILE *outputfile;
    outputfile = fopen(outputName, "w");
    fprintf(outputfile, "epsilon = %lf\n", EPSILON);
    fprintf(outputfile, "kappa = %lf\n", KAPPA);
    fprintf(outputfile, "N = %d\n", N);
    fprintf(outputfile, "NP = %d\n", NP);
    fprintf(outputfile, "dt = %lf\n", DT);
    fprintf(outputfile, "tmax = %d\n", TMAX);
    fprintf(outputfile, "NP = %d\n", NP);
    fprintf(outputfile, "L = %lf\n", L);
    fprintf(outputfile, "KMAX = %d\n", KMAX);
    fprintf(outputfile, "flow rate = %lf\n", flowRate);
    fprintf(outputfile, "\n\n");
    fprintf(outputfile, "SEED %ld\n", *idum);
    fclose(outputfile);
}

void initParaEvir(int numThreads){
    bool isOpenMP = false;
#ifdef _OPENMP
    isOpenMP = true;
    omp_set_num_threads(numThreads);
#endif
    if(isOpenMP){
        openblas_set_num_threads(1);
    } else{
        openblas_set_num_threads(numThreads);
    }
    openblas_set_num_threads(numThreads);
}

void initChain(){
    Chain = calloc(N*NP, sizeof(Chain_t));
    
    for(int i = 0; i<NP; ++i){
        int test = 0;
        while(test==0){
            test = 1;
            Chain[i*N].rx = ran1(idum)*L;
            Chain[i*N].ry = ran1(idum)*L;
            Chain[i*N].rz = ran1(idum)*L;
            
            if(usePBC){
                Chain[i*N].rx -= round(Chain[i*N].rx/L)*L;
                Chain[i*N].ry -= round(Chain[i*N].ry/L)*L;
                Chain[i*N].rz -= round(Chain[i*N].rz/L)*L;
            }
            
            if(CheckOverlap){
                for(int j = 0; j<i; ++j){
                    double dx = Chain[i*N].rx-Chain[j*N].rx;
                    double dy = Chain[i*N].ry-Chain[j*N].ry;
                    double dz = Chain[i*N].rz-Chain[j*N].rz;
                    if(usePBC){
                        dx -= round(dx/L)*L;
                        dy -= round(dy/L)*L;
                        dz -= round(dz/L)*L;
                    }
                    if(dx*dx+dy*dy+dz*dz<4.0){
                        test = 0;
                    }
                }
            }
        }
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            int test = 0;
            while(test==0){
                test = 1;
                double theta = ran1(idum)*2.0*3.14158;
                double phi = acos(2.0*ran1(idum)-1.0);
                Chain[i*N+j].rx = Chain[i*N+j-1].rx + 2.05*cos(theta)*sin(phi);
                Chain[i*N+j].ry = Chain[i*N+j-1].ry + 2.05*sin(theta)*sin(phi);
                Chain[i*N+j].rz = Chain[i*N+j-1].rz + 2.05*cos(phi);
                
                if(usePBC){
                    Chain[i*N+j].rx -= round(Chain[i*N+j].rx/L)*L;
                    Chain[i*N+j].ry -= round(Chain[i*N+j].ry/L)*L;
                    Chain[i*N+j].rz -= round(Chain[i*N+j].rz/L)*L;
                }
                
                if(CheckOverlap){
                    for(int k = 0; k<i*N+j; ++k){
                        double dx = Chain[i*N+j].rx-Chain[k].rx;
                        double dy = Chain[i*N+j].ry-Chain[k].ry;
                        double dz = Chain[i*N+j].rz-Chain[k].rz;
                        if(usePBC){
                            dx -= round(dx/L)*L;
                            dy -= round(dy/L)*L;
                            dz -= round(dz/L)*L;
                        }
                        if(dx*dx+dy*dy+dz*dz<4.0){
                            test = 0;
                        }
                    }
                }
            }
        }
    }

}

void initRingStructure(){
    int test = 0;
    double instKappa = 0.5;
    
    while(test==0){
        for(int i = 0; i<N; ++i){
            int j;
            if(i==0) j = N-1;
            else j = i-1;
            
            double dx = Chain[i].rx - Chain[j].rx;
            double dy = Chain[i].ry - Chain[j].ry;
            double dz = Chain[i].rz - Chain[j].rz;
            
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double Fs = -instKappa*(r-2.0);
            
            Chain[i].fx += Fs*dx/r;
            Chain[i].fy += Fs*dy/r;
            Chain[i].fz += Fs*dz/r;
            Chain[j].fx -= Fs*dx/r;
            Chain[j].fy -= Fs*dy/r;
            Chain[j].fz -= Fs*dz/r;
            
        }
        
        for(int i = 0; i<N; ++i){
            Chain[i].rx += DT*Chain[i].fx;
            Chain[i].ry += DT*Chain[i].fy;
            Chain[i].rz += DT*Chain[i].fz;
        }
        
        double dx = Chain[N-1].rx - Chain[0].rx;
        double dy = Chain[N-1].ry - Chain[0].ry;
        double dz = Chain[N-1].rz - Chain[0].rz;
        
        if(dx*dx+dy*dy+dz*dz <= 5.0){
            test = 1;
        }
    }
}

void initM2(){
    int m2d = KMAX*2.0 + 1.0;
    m2 = calloc(m2d, sizeof(m2[0]));
    for(int i = 0; i < m2d; i++){
        m2[i] = calloc(m2d, sizeof(m2[0][0]));
        for(int j = 0; j < m2d; j++){
            m2[i][j] = calloc(m2d, sizeof(m2[0][0][0]));
            for(int k = 0; k < m2d; k++){
                m2[i][j][k] = calloc(9, sizeof(m2[0][0][0][0]));
            }
        }
    }
    
    double alpha = 6.0/L;
    
//#pragma omp parallel for schedule (static)
    for(int kx = -KMAX; kx<KMAX+1; ++kx){
        double rkx = 2.0*M_PI*kx/L;
        for(int ky = -KMAX; ky<KMAX+1; ++ky){
            double rky = 2.0*M_PI*ky/L;
            for(int kz = -KMAX; kz<KMAX+1; ++kz){
                double rkz = 2.0*M_PI*kz/L;
                double rkk = rkx*rkx+rky*rky+rkz*rkz;
                double mm3 = (1.0-rkk/3.0)*(1.0+rkk*pow(alpha,-2.0)/4.0+rkk*rkk*pow(alpha,-4.0)/8.0)*6.0*M_PI/rkk*exp(-rkk*pow(alpha,-2.0)/4.0);
                
                m2[kx+KMAX][ky+KMAX][kz+KMAX][0] += mm3*(1.0-rkx*rkx/rkk);
                m2[kx+KMAX][ky+KMAX][kz+KMAX][1] += -mm3*rkx*rky/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][2] += -mm3*rkx*rkz/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][3] += -mm3*rky*rkx/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][4] += mm3*(1.0-rky*rky/rkk);
                m2[kx+KMAX][ky+KMAX][kz+KMAX][5] += -mm3*rky*rkz/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][6] += -mm3*rkz*rkx/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][7] += -mm3*rkz*rky/rkk;
                m2[kx+KMAX][ky+KMAX][kz+KMAX][8] += mm3*(1.0-rkz*rkz/rkk);
            }
        }
    }
}

long initRan(){
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000+1000000000;
}

float ran1(long *idum){
    long j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if(*idum <= 0){
        if(-(*idum)<1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void printInitial(){
    FILE *Trajectory;
    Trajectory = fopen(trajName,"a");
    fprintf(Trajectory, "%d\n-1\n", N*NP);
    for(int i = 0; i<N*NP; ++i){
        fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}


