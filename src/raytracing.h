/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   raytracing.h
 * Author: davidjwilliamson
 *
 * Created on 19 April 2021, 14:29
 */

#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <memory>
#include "coupling.h"

struct Pos_Type {
    double Pos[3];
};

class Raytracer {
    public:
        Raytracer(std::shared_ptr<ParticlePositionCoupler> particles);
        
        void build_tree();
        void agn_optical_depths(double* r_agn, double *depths, bool agn_at_origin);

    private:

        double one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg, int localOffset,int localNActive,int localIndex[]
                                ,double *tdiffbeam);
        double one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg, int localOffset,int localNActive,int localIndex[]
                                ,double *tdiffbeam, double r_agn[]);
        double one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg, int localOffset,int localNActive,int localIndex[]
                                ,double *tdiffbeam, double r_agn[],bool agn_at_origin);

        int packParticles(struct Pos_Type MyPos[], int localIndex[]);
        int gatherParticles(struct Pos_Type AllPos[], struct Pos_Type MyPos[], int numpartlist[], int numpartdisplacements[], int localNActive, double *tdiffwait, double *tdiffcomm);

        void setup_agn_kernel();
        
        std::shared_ptr<AGN_Kernel> agn_kernel = nullptr;
        std::shared_ptr<absorbtree> tree = nullptr;
        std::shared_ptr<ParticlePositionCoupler> particles = nullptr;

        std::vector<int>** lol;
        const unsigned int maxlolsize = 1000000;
    
};

void setup_raytracing();



#endif /* RAYTRACING_H */

