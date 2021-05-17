/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ArrayParticlePositionCoupler.cpp
 * Author: davidjwilliamson
 * 
 * Created on 28 April 2021, 14:51
 */

#include "ArrayParticlePositionCoupler.h"


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
