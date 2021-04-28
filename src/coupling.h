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


#ifndef COUPLING_H
#define COUPLING_H

#include <vector>

class ParticlePositionCoupler {
public:
    virtual double *position(int i) {};
    virtual double mass(int i) {};
    virtual double smoothing(int i) {};
    virtual double opacity(int i) {};
    virtual bool isDusty(int i) {};
    virtual bool isGas(int i) {};
    virtual bool isActive(int i) {};
    virtual bool isAlive(int i) {};
    virtual long long allN() {};
    virtual long long localN() {};
    virtual long long allGasN() {};
};

//extern ArrayParticlePositionCoupler *particlePositionCoupler;

double opticalDepthCutoff();

double kernel_wk(double r,double h);

double dtime_system(double t0, double t1);

double system_time(void);

//    void endrun(int code) {
//        
////        std::cout << "dying, error " << code << std::endl;
////        exit(0); // do MPI finalize first?
//    }

#endif /* COUPLING_H */

