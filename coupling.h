/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   coupling.h
 * Author: davidjwilliamson
 *
 * Created on 13 April 2021, 16:25
 */
//#ifdef __cplusplus
//extern "C" {
//#endif
//#ifdef __cplusplus
//}
//#endif

#include "../../allvars.h"
#include "../../kernel.h"

extern "C"
{ 
#include "../../proto.h"
}

#ifndef COUPLING_H
#define COUPLING_H

#ifdef __cplusplus
extern "C" {
#endif


    double inline *position(int i) {
        return P[i].Pos;
    }

    double inline mass(int i) {
        return P[i].Mass;
    }
    
    double inline smoothing(int i) {
        return P[i].Hsml;
    }

    double inline opacity(int i) {
        return SphP[i].Opacity;
    }
    
    bool inline isDusty(int i) {
        // non-sputtered criterion
        return SphP[i].InternalEnergy<=All.u_DustSput;
    }
    
    bool inline isGas(int i) {
        return P[i].Type==0;
    }
    
    bool inline isActive(int i) {
        return TimeBinActive[P[i].TimeBin];
    }

    bool inline isAlive(int i) {
        return P[i].Mass>0.;
    }
    
    double inline opticalDepthCutoff() {
        return All.AGNOpticalDepthCutoff;
    }
    
    long long inline allGasN() {
        return All.TotN_gas;
    }

    double inline kernel_wk(double r,double h) {
        if ( r>=h ) return 0.;

        // get kernel from GIZMO
        double hinv, hinv3, hinv4,wk,dwk;
        kernel_hinv(h, &hinv, &hinv3, &hinv4);
        kernel_main(r/h, hinv3, hinv4, &wk, &dwk, 0);

        return wk;        
    }
    
    
    double inline dtime_system(double t0, double t1) {
        return timediff(t0,t1);
    }
    
    double inline system_time(void) {
        return my_second();
    }

#ifdef __cplusplus
}
#endif

#endif /* COUPLING_H */

