/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ArrayParticlePositionCoupler.h
 * Author: davidjwilliamson
 *
 * Created on 28 April 2021, 14:51
 */

#ifndef ARRAYPARTICLEPOSITIONCOUPLER_H
#define ARRAYPARTICLEPOSITIONCOUPLER_H
#include "coupling.h"

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


#endif /* ARRAYPARTICLEPOSITIONCOUPLER_H */

