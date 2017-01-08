//
//  main.c
//  SingleChain
//
//  Created by Linling Miao on 6/7/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "main.h"
#include "Initialization.h"
#include "MobilityMatrix.h"
#include "GeyerWinterCalculation.h"

int main(int argc, char * argv[]) {
    
    initialization(argc, argv);
    
    for(int t = 0; t<TMAX; ++t){
        //Starttime = timer();
        
        resetForce();
        storeCoM();
        
        forceSpring();
        forceRepulsive();
        
        if(usePBC) EwaldSumRPY();
        else regularRPY();
        
        AvgMatrix(t);
        
        if(GeyerWinter){
            GeyerWinterCalculation();
            updateChain();
        } else{
            Cholesky(t);
            updateChainCH();
        }
        
        Rgyration(t);
        LongTimeDc(t);
        Rendtoend(t);
        Extension(t);
        
        if(usePBC) applyPBC();
        
        printTrajectory(t);
        
        //long end = timer() - Starttime;
        //printf("%ld\n",end);
    }
    
    return 0;
}

Vector3D_t getNID(int i, int j){
    Vector3D_t NID;
    NID.x = Chain[i].rx - Chain[j].rx;
    NID.y = Chain[i].ry - Chain[j].ry;
    NID.z = Chain[i].rz - Chain[j].rz;
    
    NID.x -= round(NID.x/L)*L;
    NID.y -= round(NID.y/L)*L;
    NID.z -= round(NID.z/L)*L;
    
    return NID;
}

void applyPBC(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].rx -= round(Chain[i].rx/L)*L;
        Chain[i].ry -= round(Chain[i].ry/L)*L;
        Chain[i].rz -= round(Chain[i].rz/L)*L;
    }
}

void resetForce(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].fx = 0.0;
        Chain[i].fy = 0.0;
        Chain[i].fz = 0.0;
    }
}

void forceSpring(){
//#pragma omp parallel for collapse(2)
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            Vector3D_t NID;
            
            if(usePBC){
                NID = getNID(i*N+j,i*N+j-1);
            } else{
                NID.x = Chain[i*N+j].rx - Chain[i*N+j-1].rx;
                NID.y = Chain[i*N+j].ry - Chain[i*N+j-1].ry;
                NID.z = Chain[i*N+j].rz - Chain[i*N+j-1].rz;
            }
            
            double r = sqrt(NID.x*NID.x+NID.y*NID.y+NID.z*NID.z);
            double Fs = -KAPPA*(r-2.0);
            
            Chain[i*N+j].fx += Fs*NID.x/r;
            Chain[i*N+j].fy += Fs*NID.y/r;
            Chain[i*N+j].fz += Fs*NID.z/r;
            Chain[i*N+j-1].fx -= Fs*NID.x/r;
            Chain[i*N+j-1].fy -= Fs*NID.y/r;
            Chain[i*N+j-1].fz -= Fs*NID.z/r;
            
        }
    }
}

void forceRepulsive(){
//#pragma omp parallel for collapse(2)
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<N; ++j){
            for(int k = i*N+j+1; k<N*NP; ++k){
                Vector3D_t NID;
                
                if(usePBC){
                    NID = getNID(i*N+j,k);
                } else{
                    NID.x = Chain[i*N+j].rx - Chain[k].rx;
                    NID.y = Chain[i*N+j].ry - Chain[k].ry;
                    NID.z = Chain[i*N+j].rz - Chain[k].rz;
                }
                
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                
                if(rr<25){
                    double ratio = 4.00/rr;
                    double r6 = ratio*ratio*ratio;
                    if(r6>3) r6 = 3;
                    double coeff = (12.0*EPSILON/rr)*(r6*r6-r6);
                    Chain[i*N+j].fx += coeff*NID.x;
                    Chain[i*N+j].fy += coeff*NID.y;
                    Chain[i*N+j].fz += coeff*NID.z;
                    Chain[k].fx -= coeff*NID.x;
                    Chain[k].fy -= coeff*NID.y;
                    Chain[k].fz -= coeff*NID.z;
                    
                }
            }
        }
    }
}

void forceSpringRing(){
    //#pragma omp parallel for private(i) schedule(dynamic)
    for(int i = 0; i<N; ++i){
        int j;
        if(i==0) j = N-1;
        else j = i-1;
        
        Vector3D_t NID;
        
        NID.x = Chain[i].rx - Chain[j].rx;
        NID.y = Chain[i].ry - Chain[j].ry;
        NID.z = Chain[i].rz - Chain[j].rz;
        
        double r = sqrt(NID.x*NID.x+NID.y*NID.y+NID.z*NID.z);
        double Fs = -KAPPA*(r-2.0);
        
        Chain[i].fx += Fs*NID.x/r;
        Chain[i].fy += Fs*NID.y/r;
        Chain[i].fz += Fs*NID.z/r;
        Chain[j].fx -= Fs*NID.x/r;
        Chain[j].fy -= Fs*NID.y/r;
        Chain[j].fz -= Fs*NID.z/r;
    }
}

void updateChainCH(){
    int nn =  3*N;
    double p = sqrt(2.0*DT);
    
    float *RR = calloc(nn, sizeof(float));
    for(int i = 0; i<nn; ++i){
        RR[i] = gasdev(idum);
    }
    
    float *D_noise = calloc(nn,sizeof(float));
    float *D_force = calloc(nn,sizeof(float));
    float *force = calloc(nn,sizeof(float));
    
    for(int i = 0; i<N; ++i){
        force[i*3] = Chain[i].fx;
        force[i*3+1] = Chain[i].fy;
        force[i*3+2] = Chain[i].fz;
    }
    
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,matrix->mm,nn,force,1,1,D_force,1);
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,lm,nn,RR,1,1,D_noise,1);
    
    /*Update Position*/
    if(flowType == 0){ //No flow
        for(int i = 0; i<N; ++i){
            Chain[i].rx += DT*D_force[i*3] + p*D_noise[i*3];
            Chain[i].ry += DT*D_force[i*3+1] + p*D_noise[i*3+1];
            Chain[i].rz += DT*D_force[i*3+2] + p*D_noise[i*3+2];
        }
    } else if(flowType == 1){ //elongational flow
        for(int i = 0; i<N; ++i){
            Chain[i].rx += DT*D_force[i*3] + DT*flowRate*(Chain[i].rx - R_CoM[0].x) + p*D_noise[i*3];
            Chain[i].ry += DT*D_force[i*3+1] - DT*flowRate*(Chain[i].ry - R_CoM[0].y) + p*D_noise[i*3+1];
            Chain[i].rz += DT*D_force[i*3+2] + p*D_noise[i*3+2];
        }
    } else if(flowType == 2){ //steady shear
        for(int i = 0; i<N; ++i){
            Chain[i].rx += DT*D_force[i*3] + DT*flowRate*(Chain[i].ry - R_CoM[0].y) + p*D_noise[i*3];
            Chain[i].ry += DT*D_force[i*3+1] + p*D_noise[i*3+1];
            Chain[i].rz += DT*D_force[i*3+2] + p*D_noise[i*3+2];
        }
    }
    
    free(force);
    free(D_force);
    free(D_noise);
    free(RR);
}

void updateChain(){
    int nn =  3*N;
    double p = sqrt(2.0*DT);
    
    float *RR = calloc(nn, sizeof(float));
    
//#pragma omp parallel for chedule(static)
    for(int i = 0; i<NP; ++i){
        float* D_noise = calloc(nn,sizeof(float));
        float* D_force = calloc(nn,sizeof(float));
        float* D_ij = calloc(nn*nn,sizeof(float));
        float* force = calloc(nn,sizeof(float));
    
        for(int j = 0; j<NP; ++j){
            locateD_ij(i,j,D_ij);
            locateForce(j,force);
            
            for(int k = 0; k<nn; ++k){
                RR[k] = gasdev(idum);
            }
            
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,D_ij,nn,force,1,1,D_force,1);
            /*Since beta_ii = 1, devided by beta_ij before multiply it*/
            for(int l = 0; l<nn; ++l){
                D_ij[l*nn+l] /= gwParm->beta_ij;
            }
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,gwParm->beta_ij,D_ij,nn,RR,1,1,D_noise,1);
        }
        
        /*Update Position*/
        if(flowType == 0){ //No flow
            for(int j = 0; j<N; ++j){
                Chain[i*N+j].rx += DT*D_force[j*3] + p*gwParm->C_i[j*3]*D_noise[j*3];
                Chain[i*N+j].ry += DT*D_force[j*3+1] + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                Chain[i*N+j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
            }
        } else if(flowType == 1){ //elongational flow
            for(int j = 0; j<N; ++j){
                Chain[i*N+j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[i*N+j].rx - R_CoM[i].x) + p*gwParm->C_i[j*3]*D_noise[j*3];
                Chain[i*N+j].ry += DT*D_force[j*3+1] - DT*flowRate*(Chain[i*N+j].ry - R_CoM[i].y) + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                Chain[i*N+j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
            }
        } else if(flowType == 2){ //steady shear
            for(int j = 0; j<N; ++j){
                Chain[i*N+j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[i*N+j].ry - R_CoM[i].y) + p*gwParm->C_i[j*3]*D_noise[j*3];
                Chain[i*N+j].ry += DT*D_force[j*3+1] + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                Chain[i*N+j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
            }
        }
        
        
        free(force);
        free(D_ij);
        
        free(D_noise);
        free(D_force);
        
    }
    
    free(RR);
}

void locateD_ij(int i, int j, float *D_ij){
    int nn = 3*N;
    int start = i*NP*nn*nn + j*nn*nn;
    for(int k = 0; k<nn*nn; ++k){
        D_ij[k] = matrix->mm[start+k];
    }
}

void locateForce(int j, float *force){
    for(int k = 0; k<N; ++k){
        force[k*3] = Chain[j*N+k].fx;
        force[k*3+1] = Chain[j*N+k].fy;
        force[k*3+2] = Chain[j*N+k].fz;
    }
}

Vector3D_t CenterOfMass(int i){
    Vector3D_t CoM;
    if(usePBC){
        CoM.x = Chain[i*N].rx;
        CoM.y = Chain[i*N].ry;
        CoM.z = Chain[i*N].rz;
        for(int j = 1; j<N; ++j){
            Vector3D_t NID = getNID(i*N+j,i*N+j-1);
            CoM.x += (N-j)*NID.x/N;
            CoM.y += (N-j)*NID.y/N;
            CoM.z += (N-j)*NID.z/N;
        }
        CoM.x -= round(CoM.x/L)*L;
        CoM.y -= round(CoM.y/L)*L;
        CoM.z -= round(CoM.z/L)*L;
    } else{
        CoM.x = Chain[i*N].rx/N;
        CoM.y = Chain[i*N].ry/N;
        CoM.z = Chain[i*N].rz/N;
        for(int j = 1; j<N; ++j){
            CoM.x += Chain[i*N+j].rx/N;
            CoM.y += Chain[i*N+j].ry/N;
            CoM.z += Chain[i*N+j].rz/N;
        }
    }
    return CoM;
}

void storeCoM(){
    if(R_CoM==NULL){
        R_CoM = calloc(NP,sizeof(Vector3D_t));
    }
    
    for(int i = 0; i<NP; ++i){
        R_CoM[i] = CenterOfMass(i);
    }
}

void LongTimeDc(int t){
    for(int i = 0; i<NP; ++i){
        Vector3D_t dCoM;
        Vector3D_t new_CoM = CenterOfMass(i);
        dCoM.x = new_CoM.x - R_CoM[i].x;
        dCoM.y = new_CoM.y - R_CoM[i].y;
        dCoM.z = new_CoM.z - R_CoM[i].z;
        
        if(usePBC){
            dCoM.x -= round(dCoM.x/L)*L;
            dCoM.y -= round(dCoM.y/L)*L;
            dCoM.z -= round(dCoM.z/L)*L;
        }
        
        MSD_CoM += dCoM.x*dCoM.x+dCoM.y*dCoM.y+dCoM.z*dCoM.z;
    }
}

void Rgyration(int t){
    if(t==0){
        Rg = 0.0; Ree = 0.0;
    }
    
    //Radius of Gyration
    for(int i = 0; i<NP; ++i){
        Vector3D_t CoM = CenterOfMass(i);
        //printf("%lf %lf %lf\n",CoM.x,CoM.y,CoM.z);
        for(int j = 0; j<N; ++j){
            double dx = Chain[i*N+j].rx - CoM.x;
            double dy = Chain[i*N+j].ry - CoM.y;
            double dz = Chain[i*N+j].rz - CoM.z;
            
            Rg += dx*dx+dy*dy+dz*dz;
        }
        
        /*if(t>0 && t%1000==0){
            FILE *diffusivity = fopen(DiffName,"a");
            fprintf(diffusivity, "%lf \n",MSD_CoM*1000/6.0/t/NP);
            fclose(diffusivity);
        }*/
    }
    
    //End to End Vector
    if(usePBC){
        for(int i = 0; i<NP; ++i){
            double dx = 0;
            double dy = 0;
            double dz = 0;
            for(int j = 1; j<N; ++j){
                Vector3D_t NID = getNID(i*N+j,i*N+j-1);
                dx += NID.x;
                dy += NID.y;
                dz += NID.z;
            }
            Ree += dx*dx+dy*dy+dz*dz;
        }
    } else{
        for(int i = 0; i<NP; ++i){
            double dx = Chain[i*N].rx - Chain[(i+1)*N-1].rx;
            double dy = Chain[i*N].ry - Chain[(i+1)*N-1].ry;
            double dz = Chain[i*N].rz - Chain[(i+1)*N-1].rz;
            
            Ree += dx*dx+dy*dy+dz*dz;
        }
    }
    
    /*if(t>0 && t%1000==0){
        FILE *statproperty = fopen(RgName,"a");
        fprintf(statproperty,"%lf %lf\n",Rg/1000/NP/N, Ree/1000/NP);
        fclose(statproperty);
    }*/
    
    if(t%1000==0){
        Rg = 0.0; Ree = 0.0;
    }
}

void Rendtoend(int t){
    double dx = 0;
    double dy = 0;
    double dz = 0;
    if(usePBC) {
        for(int i = 1; i<N; ++i){
            Vector3D_t NID = getNID(i,i-1);
            dx += NID.x;
            dy += NID.y;
            dz += NID.z;
        }
    } else{
        dx = Chain[N-1].rx - Chain[0].rx;
        dy = Chain[N-1].ry - Chain[0].ry;
        dz = Chain[N-1].rz - Chain[0].rz;
    }
    
    /*if(t>0 && t%1000==0){
        FILE *ree = fopen(reeName,"a");
        fprintf(ree,"%lf %lf %lf\n", dx,dy,dz);
        fclose(ree);
    }*/
}
void Extension(int t){
    for(int i = 0; i<NP; ++i){
        double dxmax = 0.0;
        double dxmin = 0.0;
        double ddx = 0.0;
        
        for(int j = 1; j<N; ++j){
            Vector3D_t NID = getNID(i*N+j,i*N+j-1);
            ddx += NID.x;
            if(ddx<dxmin) dxmin = ddx;
            else if(ddx>dxmax) dxmax = ddx;
        }
        ext.x += dxmax - dxmin;
    }
    
    /*if(t>0 && t%1000==0){
        FILE *Extensionout;
        Extensionout = fopen(extName, "a");
        fprintf(Extensionout, "%lf\n",ext/1000/NP);
        fclose(Extensionout);
    }*/
    
    if(t%1000==0) ext.x = 0.0;
}
void ExtensionRing(int t){
    for(int i = 0; i<NP; ++i){
        double dxmax = 0.0; double dymax = 0.0; double dzmax = 0.0;
        double dxmin = 0.0; double dymin = 0.0; double dzmin = 0.0;
        double ddx = 0.0; double ddy = 0.0; double ddz = 0.0;
        
        for(int j = 1; j<N; ++j){
            Vector3D_t NID;
            NID.x = Chain[i*N+j].rx - Chain[i*N+j-1].rx;
            NID.y = Chain[i*N+j].ry - Chain[i*N+j-1].ry;
            NID.z = Chain[i*N+j].rz - Chain[i*N+j-1].rz;
            
            ddx += NID.x;
            ddy += NID.y;
            ddz += NID.z;
            if(ddx<dxmin) dxmin = ddx;
            else if(ddx>dxmax) dxmax = ddx;
            if(ddy<dymin) dymin = ddy;
            else if(ddy>dymax) dymax = ddy;
            if(ddz<dzmin) dzmin = ddz;
            else if(ddz>dzmax) dzmax = ddz;
        }
        ext.x += dxmax - dxmin;
        ext.y += dymax - dymin;
        ext.z += dzmax - dzmin;
    }
    
    /*if(t>0 && t%1000==0){
        FILE *Extensionout;
        Extensionout = fopen(extName, "a");
        fprintf(Extensionout, "%lf %lf %lf\n", ext.dx/1000, ext.dy/1000, ext.dz/1000);
        fclose(Extensionout);
    }*/
    
    if(t%1000==0){
        ext.x = 0.0; ext.y = 0.0; ext.z = 0.0;
    }
}

float gasdev(long *idum){
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0){
        do{
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else{
        iset = 0;
        return gset;
    }
}

void printTrajectory(int t){
    FILE *Trajectory;
    if(t%trajStep==0){
        Trajectory = fopen(trajName,"a");
        fprintf(Trajectory, "%d\n%d\n", N*NP, t);
        int i;
        for(i = 0; i<N*NP; ++i){
            fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
        }
        fclose(Trajectory);
    }
}



long long timer(){
    struct timespec ts;
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return 1000000000*ts.tv_sec + ts.tv_nsec;
}
