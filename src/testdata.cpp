/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testdata.cpp
 * Author: davidjwilliamson
 * 
 * Created on 16 April 2021, 10:09
 */

#include "coupling.h"
#include "testdata.h"
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <random>
#include <memory>

// globals, from testdata.h
//struct TestData *testpositions;
//int Nlocal,Ntot;
//ArrayParticlePositionCoupler *particlePositionCoupler = new ArrayParticlePositionCoupler();

std::shared_ptr<ArrayParticlePositionCoupler> generate_test_data(int N) {
    
    std::shared_ptr<ArrayParticlePositionCoupler> particlePositionCoupler(new ArrayParticlePositionCoupler());
    
    int NTask,ThisTask;
    
    MPI_Comm_size(MPI_COMM_WORLD, &NTask); 
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); 
    
    particlePositionCoupler->Nlocal = N;
    particlePositionCoupler->Ntot = N * NTask;

    particlePositionCoupler->testpositions = new TestData[N];
    
    // creates a cube, NTask**(1/3) per side
    // if NTask isn't a cube, there will be some spillover
    // whatever, this is just a quick test
    
    int nside = std::lround(std::cbrt(NTask));
    
    int offset[3] = {ThisTask%nside,(ThisTask/nside)%nside,ThisTask/(nside*nside)};
    
    double cent = -((float)nside)/2.;
    
//    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(ThisTask); //Standard mersenne_twister_engine seeded with proc number for consistency
    std::normal_distribution<> dis(0.5, 0.2);
    
    for ( int ip=0 ; ip<N ; ip++ ) {
        for ( int jj=0 ; jj<3 ; jj++ ) {
            particlePositionCoupler->testpositions[ip].Pos[jj] = dis(gen)+offset[jj]+cent;
        }
        // flatten for easier viz, harder tree
        particlePositionCoupler->testpositions[ip].Pos[2]*=0.01;
    }
    
    return particlePositionCoupler;
    
//    std::cout << ThisTask << " " << ix << " " << iy << " " << iz << " " << nside << std::endl;
}

void dump_positions(std::shared_ptr<ArrayParticlePositionCoupler> particlePositionCoupler) {
    int NTask,ThisTask;
    
    MPI_Comm_size(MPI_COMM_WORLD, &NTask); 
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); 
    
    TestData *allData;
    
    if (ThisTask==0) {
        allData = new TestData[particlePositionCoupler->allN()];
    }
    
    MPI_Gather(particlePositionCoupler->testpositions, 
           sizeof(TestData)*particlePositionCoupler->localN(), 
           MPI_BYTE,
           allData, 
           sizeof(TestData)*particlePositionCoupler->localN(), 
           MPI_BYTE, 
           0, 
           MPI_COMM_WORLD);
    
    if (ThisTask==0) {
        for ( int ip=0 ; ip<particlePositionCoupler->allN() ; ip++ ) {
            for ( int jj=0 ; jj<3 ; jj++ ) {
                std::cout << allData[ip].Pos[jj] <<" ";
            }
            std::cout << allData[ip].OpticalDepth << std::endl << std::flush;
        }
        delete[] allData;
    }
}

