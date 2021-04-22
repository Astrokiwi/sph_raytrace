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

struct TestData {
    double Pos[3];
    double OpticalDepth;
    double smoothing;
    double opacity;
    double mass;
    int proc;
};

class ArrayParticlePositionCoupler: public ParticlePositionCoupler {
public:
    std::vector<struct TestData> testpositions;

    int Nlocal,Ntot;
    
    double constant_mass = 1.;

    double *position(int i) override;
    double mass(int i) override;
    double smoothing(int i) override;
    double opacity(int i) override;
    bool isDusty(int i) override;
    bool isGas(int i) override;
    bool isActive(int i) override;
    bool isAlive(int i) override;
    long long allN() override;
    long long localN() override;
    long long allGasN() override;
    
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

