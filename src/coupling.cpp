/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "coupling.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "testdata.h"

double *ArrayParticlePositionCoupler::position(int i) {
    return this->testpositions[i].Pos;
}

double ArrayParticlePositionCoupler::mass(int i) {
    return constant_mass;
}

double ArrayParticlePositionCoupler::smoothing(int i) {
    return this->testpositions[i].smoothing;
}

double ArrayParticlePositionCoupler::opacity(int i) {
    return this->testpositions[i].opacity;
}

bool ArrayParticlePositionCoupler::isDusty(int i) {
    // non-sputtered criterion
    return true;
}

bool ArrayParticlePositionCoupler::isGas(int i) {
    return true;
}

bool ArrayParticlePositionCoupler::isActive(int i) {
    return true;
}

bool ArrayParticlePositionCoupler::isAlive(int i) {
    return true;
}

long long ArrayParticlePositionCoupler::allN() {
    return this->Ntot;
}

long long ArrayParticlePositionCoupler::localN() {
    return this->Nlocal;
}

long long ArrayParticlePositionCoupler::allGasN() {
    return this->Ntot;
}


double opticalDepthCutoff() {
    return 7e6;
}
double kernel_wk(double r,double h) {
    if ( r>=h ) return 0.;



    double hinv, hinv3, wk;

    double u = r/h;

    // kernel ripped & butchered from GIZMO
    hinv = 1.0 / h;

    hinv3 = hinv * hinv * hinv;

    if(u < 0.5)
      {
            wk = (1.0 + 6.0 * (u - 1.0) * u * u);
      }
    else
      {
        double t1 = (1.0 - u);
        double t2 = t1 * t1;
            wk = 2.0 * t2 * t1;
      }



    wk *= (15.0/(2.0*M_PI)) * hinv3;

    return wk;        
}


double dtime_system(double t0, double t1) {
    // TODO: implement
    return 0.;
}

double system_time(void) {
    // TODO: implement
    return 0.;
}
