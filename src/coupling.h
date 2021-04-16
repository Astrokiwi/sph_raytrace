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

double *position(int i);

double mass(int i);

double smoothing(int i);

double opacity(int i);

bool isDusty(int i);

bool isGas(int i);

bool isActive(int i);

bool isAlive(int i);

double opticalDepthCutoff();

long long allN();

long long localN();

long long allGasN();

double kernel_wk(double r,double h);

double dtime_system(double t0, double t1);

double system_time(void);

//    void endrun(int code) {
//        
////        std::cout << "dying, error " << code << std::endl;
////        exit(0); // do MPI finalize first?
//    }

#endif /* COUPLING_H */

