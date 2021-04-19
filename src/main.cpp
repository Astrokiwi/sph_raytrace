/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: davidjwilliamson
 *
 * Created on 15 April 2021, 13:48
 */

#include <cstdlib>
#include <mpi.h>
#include "rt_prototypes.h"
#include "testdata.h"
#include "coupling.h"
#include "raytracing.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    
    std::shared_ptr<ArrayParticlePositionCoupler> p = generate_test_data(200);
    
    std::unique_ptr<Raytracer> raytracer(new Raytracer(p));
    raytracer->build_tree();

    double *AGN_localtau = new double[p->localN()];
    double r_agn[3] = {0,0,0};
    raytracer->agn_optical_depths(r_agn,AGN_localtau,true);

    for ( int ii=0 ; ii<p->localN() ; ii++ ) {
        p->testpositions[ii].OpticalDepth = AGN_localtau[ii];
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    dump_positions(p);

    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();		/* clean up & finalize MPI */

    return 0;
}

