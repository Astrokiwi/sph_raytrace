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
#include "rt_prototypes.h"
#include "raytracing.h"
#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <random>
#include <memory>

// globals, from testdata.h
//struct TestData *testpositions;
//int Nlocal,Ntot;
//ArrayParticlePositionCoupler *particlePositionCoupler = new ArrayParticlePositionCoupler();

std::shared_ptr<ArrayParticlePositionCoupler> distribute_data(std::vector<TestData> &allData, int Ntot) {
    // Distribute global data
    
    int NTask,ThisTask;
    
    std::shared_ptr<ArrayParticlePositionCoupler> particlePositionCoupler(new ArrayParticlePositionCoupler());
    particlePositionCoupler->Ntot = Ntot;

    
    MPI_Comm_size(MPI_COMM_WORLD, &NTask); 
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); 

    {
        int counts[NTask],displ[NTask];
        
        int n_mean = particlePositionCoupler->Ntot/NTask;

        // Local allocations

    //    particlePositionCoupler->testpositions = new TestData[N];
        
        for ( int iTask=0 ; iTask<NTask-1 ; iTask++ ) {
            counts[iTask]=n_mean;
        }
        counts[NTask-1]=particlePositionCoupler->Ntot-n_mean*(NTask-1); // last one gets the remainder

        particlePositionCoupler->Nlocal = counts[ThisTask];
        particlePositionCoupler->testpositions.resize(particlePositionCoupler->Nlocal);


        for ( int iTask=0 ; iTask<NTask ; iTask++ ) {
//            if ( ThisTask==0 ) {
//                std::cout << iTask << " " << counts[iTask] << std::endl;
//            }
            counts[iTask]*=sizeof(TestData);
        }

        displ[0] = 0;
        for ( int iTask=1 ; iTask<NTask ; iTask++ ) {
            displ[iTask]=displ[iTask-1]+counts[iTask-1];
        }

        // Distribute between tasks
        MPI_Scatterv(allData.data(), 
               counts,
               displ,
               MPI_BYTE,
               particlePositionCoupler->testpositions.data(), 
               sizeof(TestData)*particlePositionCoupler->localN(), 
               MPI_BYTE, 
               0, 
               MPI_COMM_WORLD);
    }
    
    for ( int ip=0 ; ip<particlePositionCoupler->localN() ; ip++ ) {
        particlePositionCoupler->testpositions[ip].proc = ThisTask;
    }    


    return particlePositionCoupler;
}

std::shared_ptr<ArrayParticlePositionCoupler> generate_test_data() {
    
    int nblob = 27;
    int np_blob = 200;
//    int nblob = 2;
//    int np_blob = 32;
    int ThisTask;

    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); 


    // Create global data

    // creates a cube, NTask**(1/3) per side
    // if NTask isn't a cube, there will be some spillover
    // whatever, this is just a quick test
    std::vector<TestData> allData;

    if ( ThisTask==0 ) {
        int nside = std::lround(std::cbrt(nblob));
        double cent = -((float)nside)/2.;

        allData.resize(np_blob * nblob);
//        allData = new TestData[particlePositionCoupler->Ntot];

    //    std::random_device rd;  //Will be used to obtain a seed for the random number engine

        for ( int iblob=0 ; iblob<nblob ; iblob++ ) {
            int offset[3] = {iblob%nside,(iblob/nside)%nside,iblob/(nside*nside)};

            std::mt19937 gen(iblob); //Standard mersenne_twister_engine seeded with proc number for consistency
            std::normal_distribution<> dis(0.5, 0.2);

            for ( int ip=np_blob*iblob ; ip<np_blob*(iblob+1) ; ip++ ) {
                allData[ip].smoothing = 0.1;
                allData[ip].opacity = 1.e-6;
                for ( int jj=0 ; jj<3 ; jj++ ) {
                    allData[ip].Pos[jj] = dis(gen)+offset[jj]+cent;
                }
                // flatten for easier viz, harder tree
                allData[ip].Pos[2]*=0.01;
            }
        }
    }

    
    
//    if ( ThisTask==0 ) {
//        delete[] allData;
//    }

    
    return distribute_data(allData,np_blob * nblob);
    
//    std::cout << ThisTask << " " << ix << " " << iy << " " << iz << " " << nside << std::endl;
}

void dump_positions(std::shared_ptr<ArrayParticlePositionCoupler> particlePositionCoupler) {
    int NTask,ThisTask;
    
    MPI_Comm_size(MPI_COMM_WORLD, &NTask); 
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); 
    
    std::vector<TestData> allData;
    
    if (ThisTask==0) {
        allData.resize(particlePositionCoupler->allN());
    }
    
    
    {
        int counts[NTask],displ[NTask];
//        if ( ThisTask==0 ) {
//            counts = new int[NTask];
//            displ = new int[NTask];
//        }
        int n = particlePositionCoupler->localN();
//        MPI_Barrier(MPI_COMM_WORLD);
//        std::cout << n << " " << ThisTask << std::endl << std::flush;
        MPI_Gather(&n,
                1,
                MPI_INT,
                counts,
                1,
                MPI_INT,
                0,
                MPI_COMM_WORLD
                );
//        MPI_Barrier(MPI_COMM_WORLD);
        if ( ThisTask==0 ) {
            for ( int iTask=0 ; iTask<NTask ; iTask++ ) {
                counts[iTask]*=sizeof(TestData);
            }

            displ[0] = 0;
            for ( int iTask=1 ; iTask<NTask ; iTask++ ) {
                displ[iTask]=displ[iTask-1]+counts[iTask-1];
            }

//            for ( int iTask=0 ; iTask<NTask ; iTask++ ) {
//                std::cout << counts[iTask] << " " << displ[iTask] << std::endl;
//            }
//            std::cout << std::flush;

        }
//        MPI_Barrier(MPI_COMM_WORLD);
          
                
                
        MPI_Gatherv(particlePositionCoupler->testpositions.data(), 
           sizeof(TestData)*particlePositionCoupler->localN(), 
           MPI_BYTE,
           allData.data(), 
           counts,
           displ,
           MPI_BYTE, 
           0, 
           MPI_COMM_WORLD);
    }
    if (ThisTask==0) {
        for ( int ip=0 ; ip<particlePositionCoupler->allN() ; ip++ ) {
            for ( int jj=0 ; jj<3 ; jj++ ) {
                std::cout << allData[ip].Pos[jj] <<" ";
            }
            std::cout << allData[ip].OpticalDepth << std::endl << std::flush;
//            std::cout << allData[ip].OpticalDepth << " " << allData[ip].proc << std::endl << std::flush;
        }
    }
}

void run_test() {
    int argc=0;
    char **argv;
    MPI_Init(&argc,&argv);
    
    std::shared_ptr<ArrayParticlePositionCoupler> p = generate_test_data();
    
//    std::unique_ptr<Raytracer> raytracer(new Raytracer(p));
    Raytracer raytracer(p);
    raytracer.build_tree();

    double *AGN_localtau = new double[p->localN()];
    double r_agn[3] = {0,0,0};
    
    raytracer.agn_optical_depths(r_agn,AGN_localtau,true);

    for ( int ii=0 ; ii<p->localN() ; ii++ ) {
        p->testpositions[ii].OpticalDepth = AGN_localtau[ii];
    }
    
    dump_positions(p);

    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();		/* clean up & finalize MPI */

}