#include <mpi.h>
#include <stdlib.h>
#include "coupling.h"
#include "intrinsic_depth.h"

//#include "../../allvars.h"
//extern "C"
//{ 
//#include "../../proto.h" // this is solely for dtime_system ?
//}

#include <stdio.h>
#include <limits>
#include <iomanip>

#include <iostream>
#include "rt_prototypes.h"
#include "raytrace_helpers.h"

#include <cstring>

// calc heating to this optical depth
#define HEAT_DEPTH 0.5



// fully global objects (with extern in cpp_prototypes.h)

std::vector<int>** lol;
const unsigned int maxlolsize = 1000000;
//const int maxlolsize = 1000000;

// global objects for this file

AGN_Kernel *agn_kernel;

#ifdef BENCHMARK
int n_interactions;
#endif

// struct types for mpi etc

struct Pos_Type {
    double Pos[3];
};


// prototypes for *local* functions
void setup_agn_kernel();

/*! setup AGN kernel
*/
void setup_agn_kernel() {
    agn_kernel = new AGN_Kernel();
    agn_kernel->calc_depth_optical_table(HEAT_DEPTH,mass(0),-9.,0.,90); // assumption here is that mass is constant
}

void setup_raytracing() {
    lol = new std::vector<int>*[maxlolsize];

    setup_agn_kernel();
}


// std::shared_ptr<absorbtree>

/* Builds tree, based on global distribution of particles.
*/
void Raytracer::build_tree() {

    int i,ic;

    // boundaries
    double rmin[3],rmax[3];
    double size;
    
    // calc boundaries
    for ( ic=0 ; ic<3 ; ic++ ) {
        rmin[ic] = 0.;
        rmax[ic] = 0.;
    }



    for ( i=0 ; i<NumPart ; i++ ) {
        //std::cout << ThisTask << " " << P[i].ID << " " << position(i)[2] << " " << P[i].NumNgb << std::endl;
//         position(i)[2]*=1.e-5;
//         smoothing(i) = 2.e-5;
        if( isGas(i) && isAlive(i) ) {
            for ( ic=0 ; ic<3 ; ic++ ) {
                if ( position(i)[ic]<rmin[ic] ) {
                    rmin[ic] = position(i)[ic];
                }
                if ( position(i)[ic]>rmax[ic] ) {
                    rmax[ic] = position(i)[ic];
                }
            }
        }
    }
    
    size=rmax[0]-rmin[0];
    for ( int ic=1 ; ic<3 ; ic++ ) {
        size = std::max(size,rmax[ic]-rmin[ic]);
    }
    // box contains [rmin,rmin+size)
    // so must increase size by a small amount so that the "rightmost"
    // particle actually fits strictly inside
    size*=(1.+1.e-6);

    // build tree of lists of particles for hitting smoothing lengths
    tree.reset(new absorbtree(rmin,size));
    for(i = 0; i < NumPart; i++) {
        if( isGas(i) && isAlive(i) && isDusty(i) ) {
            tree->head_addP(i,position(i),smoothing(i));
        }
    }
}

int Raytracer::packParticles(struct Pos_Type MyPos[], int localIndex[]) {
    int i,ic;

    int localNActive = 0;
//     int in_centre = 0;
    for ( i=0 ; i<NumPart ; i++ ) {
        //if( isGas(i) && isAlive(i) && TimeBinActive[P[i].TimeBin] && isDusty(i)) {

#ifdef SOTON_AGN_DOUBLESTART
        if( isGas(i) && isAlive(i) && ( isActive(i) || All.NumCurrentTiStep<=1)) {
#else
        if( isGas(i) && isAlive(i) && isActive(i) ) {
#endif

//         if( isGas(i) && isAlive(i) && TimeBinActive[P[i].TimeBin] && P[i].ID==2 ) { // TESTING ONE PARTICLE - DELETE
            // Calculate position to calculate absorption to
            // This is some optical depth into the gas particle
            // If the particle has a low enough density that it's fairly
            // optically thin, then we just calculate heating in the centre
            double skin_rad = agn_kernel->skin_rad(smoothing(i),opacity(i));
            double rad = 0.;
            for ( ic=0 ; ic<3 ; ic++ ) {
                rad+=square(position(i)[ic]);
//                 std::cout << ic << " " << position(i)[ic] << std::endl;
            }
            rad = sqrt(rad);
            

            if ( skin_rad>0. ) {

    //             std::cout << ThisTask << " " << i << " " << smoothing(i)*1.e5 << " " << skin_rad*1.e5 << " " << skin_rad/smoothing(i) << " ";
                if ( skin_rad<rad ) {
        
                    for ( ic=0 ; ic<3 ; ic++ ) {
                        MyPos[localNActive].Pos[ic]=position(i)[ic]*(1.-skin_rad/rad);
        //                 std::cout << position(i)[ic]*1.e5 << " " << position(i)[ic]*(1.-skin_rad/rad)*1.e5 << " ";
                    }
                } else {
                    for ( ic=0 ; ic<3 ; ic++ ) {
                        MyPos[localNActive].Pos[ic]=0.;
//                         in_centre++;
        //                 std::cout << position(i)[ic]*1.e5 << " " << position(i)[ic]*(1.-skin_rad/rad)*1.e5 << " ";
                    }
                }
            } else {
                for ( ic=0 ; ic<3 ; ic++ ) {
                    MyPos[localNActive].Pos[ic]=position(i)[ic];
        //                 std::cout << position(i)[ic]*1.e5 << " " << position(i)[ic]*(1.-skin_rad/rad)*1.e5 << " ";
                }
            }

//             std::cout << std::endl;
            localIndex[localNActive] = i;
            //local_skin_rad[localNActive] = skin_rad;
            localNActive++;
        }
    }
    
    return localNActive;
}

int Raytracer::gatherParticles(struct Pos_Type AllPos[], struct Pos_Type MyPos[], int numpartlist[], int numpartdisplacements[], int localNActive, double *tdiffwait, double *tdiffcomm) {
    int i;
    int maxSize;
    double tstart,tend;
    

    // Gather numbers of particles
    tstart = system_time();
    MPI_Barrier(MPI_COMM_WORLD);
    tend = system_time();
    *tdiffwait += dtime_system(tstart,tend);

    tstart = system_time();
    MPI_Allgather(&localNActive, 1, MPI_INT, &(numpartlist[0]), 1, MPI_INT, MPI_COMM_WORLD);
    tend = system_time();
    *tdiffcomm += dtime_system(tstart,tend);

    maxSize = 0;
    int nActiveTot = 0;
    for (i=0 ; i<NTask ; i++ ) {
        nActiveTot+=numpartlist[i];
        if ( numpartlist[i]>maxSize ) {
            maxSize = numpartlist[i];
        }
    }
    for (i=0 ; i<NTask ; i++ ) {
        numpartlist[i]=numpartlist[i]*sizeof(struct Pos_Type);
    }
    numpartdisplacements[0] = 0;
    for (i=1 ; i<NTask ; i++ ) {
        numpartdisplacements[i]=numpartdisplacements[i-1]+numpartlist[i-1];
    }

    tstart = system_time();
    MPI_Barrier(MPI_COMM_WORLD);
    tend = system_time();
    *tdiffwait += dtime_system(tstart,tend);

    tstart = system_time();
    MPI_Allgatherv(&(MyPos)[0], localNActive*sizeof(struct Pos_Type), MPI_BYTE, &(AllPos)[0], numpartlist, numpartdisplacements, MPI_BYTE, MPI_COMM_WORLD);
    tend = system_time();
    *tdiffcomm += dtime_system(tstart,tend);

    return nActiveTot;
}


double Raytracer::one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg,int localOffset,int localNActive,int localIndex[], double *tdiffbeam) {
    double r_agn[3] = {0,0,0};
    return one_agn_tree_ray(AllPos,alreadyDone,jg,localOffset,localNActive,localIndex,tdiffbeam,r_agn,true);
}

double Raytracer::one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg,int localOffset,int localNActive,int localIndex[], double *tdiffbeam,  double r_agn[]) {

    return one_agn_tree_ray(AllPos,alreadyDone,jg,localOffset,localNActive,localIndex,tdiffbeam,r_agn,false);
}

double Raytracer::one_agn_tree_ray(struct Pos_Type AllPos[], bool alreadyDone[],int jg,int localOffset,int localNActive,int localIndex[], double *tdiffbeam,  double r_agn[], bool agn_at_origin) {

    double tstart,tend;
    
    double *rPart_tree = &(AllPos[jg].Pos[0]);
    bool target_is_local = ( jg>=localOffset && jg<localOffset+localNActive );

    double *rAbsorb;


    unsigned int lolsize = 0;
    tstart = system_time();
    tree->beam(r_agn,rPart_tree,lol,&lolsize,ThisTask);
    tend = system_time();
    *tdiffbeam += dtime_system(tstart,tend);

    double *rPart; // adjusted for AGN centre

    if ( agn_at_origin ) {
        rPart = rPart_tree;
    } else {
        rPart = new double[3];
        for ( int ii=0 ; ii<3 ; ii++ ) {
            rPart[ii] = rPart_tree[ii]-r_agn[ii];
        }
    }

    double d12_norm2 = get_norm2(rPart);
    double d12_norm = sqrt(d12_norm2);

    memset(alreadyDone, 0, sizeof(bool)*N_gas );
    
    double depth = 0.;
//#ifdef INTRINSIC_DEPTH
    if ( target_is_local ) {
        depth = intrinsic_depth(rPart); // depths get added across processors later, but only need one proc to add small scale depth
    }
//#endif

    if ( lolsize>maxlolsize ) {
        endrun(8675309);
    }

    if ( !agn_at_origin ) {
        rAbsorb = new double[3];
    }

    for ( unsigned int ilist=0 ; ilist<lolsize ; ilist++ ) {
        std::vector<int > *hitp_set = lol[ilist];
        unsigned int hs = hitp_set->size() ;

        for ( unsigned int ii=0 ; ii<hs ; ii++ ) {
            int kg = (*hitp_set)[ii];

            if ( target_is_local ) {
                if ( localIndex[jg-localOffset]==kg ) {
                    continue;
                }
            }
            
            if ( alreadyDone[kg] ) {
                continue;
            }
            
            
            

            alreadyDone[kg] = true;
            
            if ( agn_at_origin ) {
                rAbsorb = &(position(kg)[0]);
            } else {
                for ( int ii=0 ; ii<3 ; ii++ ) {
                    rAbsorb[ii] = position(kg)[ii]-r_agn[ii];
                }
            }
                
            int bet = between_particle_agn(rPart,rAbsorb,smoothing(kg),d12_norm);
            if ( bet==0 ) continue;

            #ifdef BENCHMARK
            n_interactions++;
            #endif

            double d2int = intersect_agn_d2_nonorm(rPart,rAbsorb)/d12_norm2;
            double h2 = square(smoothing(kg));


            // within line of sight
            if ( d2int>h2 ) {
                continue;
            }
            
            if ( bet==1 ) {
                //extincting particle is fully between the end points, can use simper method
                depth+=mass(kg)*agn_kernel->flat_w2(d2int,h2)*opacity(kg);
            } else {
            // calculate line of sight
                double z3=0.;

                double source_overlap,target_overlap;

                for ( int jj=0 ; jj<3 ; jj++ ) {
                    z3+=rAbsorb[jj]*rPart[jj];
                }

                z3/=d12_norm;
            
                source_overlap = z3/smoothing(kg);
                if ( source_overlap>1. ) source_overlap=1.;
                target_overlap = (d12_norm-z3)/smoothing(kg);
                if ( target_overlap>1. ) target_overlap=1.;

                double weight_add = agn_kernel->half_flat_w2_truncated(d2int,h2,source_overlap)+agn_kernel->half_flat_w2_truncated(d2int,h2,target_overlap);
                if ( weight_add>0. ) {
                    depth+=mass(kg)*weight_add*opacity(kg);
                }
            }

            if ( depth>opticalDepthCutoff() ) {
                // depth is huge - we don't need to calculate it any further, because heating & rad pressure
                // from AGN is now negligible
                return depth;
            }
        }
    }
    if ( !agn_at_origin ) {
        delete[] rAbsorb;
        delete[] rPart;
    }
    
    return depth;
}


void Raytracer::agn_optical_depths(double* r_agn, double *depths, bool agn_at_origin) {
    // counters
    int i;
    int jp;
    
    // MPI sharing
    int numpartlist[NTask], numpartdisplacements[NTask];
    struct Pos_Type *AllPos = new struct Pos_Type[allGasN()]; // all *active* particles to be updated
    struct Pos_Type MyPos[N_gas]; // local *active* particles to have heating etc updated
    double *AGN_alltau = new double[allGasN()];
#ifdef SOTON_DUST_DUST_HEATING    
    bool alreadyDone[allGasN()];
#else
    bool alreadyDone[N_gas];
#endif

    double *local_AGN_alltau = new double[allGasN()]; // could be made smaller and more efficient
    int localOffset;
    int localIndex[N_gas];

    // timing
    double tstart,tend,t1,t0,tdiffray,tdiffcomm,tdifftree,tdiffwait,tdiffbeam;

    tdiffcomm = 0.;
    tdiffwait = 0.;
    tdiffbeam = 0.;
    t0 = system_time();

    //tree->matrix_dump(ThisTask);
    //std::cout << tree->check(10) << std::endl;

    //tree->dump();


    // Communicate points to everyone
    // pack particle positions for sending
    int localNActive = packParticles(MyPos,localIndex);
    
    int nActiveTot = gatherParticles(AllPos,MyPos,numpartlist,numpartdisplacements,localNActive,&tdiffwait,&tdiffcomm);
    
    // Calculate contribution from all of my local particles to 
    // column density from AGN to each particle using tree
    // - this is the expensive part of the routine!

    localOffset = numpartdisplacements[ThisTask]/sizeof(struct Pos_Type);

    tstart = system_time();

#ifdef BENCHMARK
    n_interactions = 0;
#endif
   if ( nActiveTot<NTask*16 ) {
       for (i=0 ; i<nActiveTot ; i++ ) {
                if ( agn_at_origin ) {
                    local_AGN_alltau[i] = one_agn_tree_ray(AllPos,alreadyDone,i,localOffset,localNActive,localIndex,&tdiffbeam);
                } else {
                    local_AGN_alltau[i] = one_agn_tree_ray(AllPos,alreadyDone,i,localOffset,localNActive,localIndex,&tdiffbeam,r_agn);
                }
        }
    } else {
        int thisChunk;
        int nChunks = std::min(8,NTask);
        
        int chunkIndices[nChunks+1];
        int chunkSize = nActiveTot/nChunks;
        for (int j=0 ; j<nChunks ; j++) {
            chunkIndices[j] = chunkSize*j;
        }
        chunkIndices[nChunks]=nActiveTot;
        
       bool    *local_cutoff_flags = new bool[nActiveTot],
               *global_cutoff_flags= new bool[nActiveTot];
       memset(local_cutoff_flags,0,nActiveTot*sizeof(bool));
       memset(global_cutoff_flags,0,nActiveTot*sizeof(bool));

       for (int j=0 ; j<nChunks ; j++) {
           thisChunk=(j+ThisTask)%nChunks;
           for (i=chunkIndices[thisChunk] ; i<chunkIndices[thisChunk+1] ; i++ ) {
               if ( !global_cutoff_flags[i] ) {
                   if ( agn_at_origin ) {
                       local_AGN_alltau[i] = one_agn_tree_ray(AllPos,alreadyDone,i,localOffset,localNActive,localIndex,&tdiffbeam);
                   } else {
                       local_AGN_alltau[i] = one_agn_tree_ray(AllPos,alreadyDone,i,localOffset,localNActive,localIndex,&tdiffbeam,r_agn);
                   }
                   if ( local_AGN_alltau[i]>opticalDepthCutoff() ) {
                    local_cutoff_flags[i] = true;
                   }
               }
           }
           MPI_Allreduce(local_AGN_alltau, AGN_alltau, nActiveTot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
           MPI_Allreduce(local_cutoff_flags, global_cutoff_flags, nActiveTot, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

//            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        delete[] local_cutoff_flags;
        delete[] global_cutoff_flags;
    }

#ifdef BENCHMARK
    {
        int globalChecksum=0;
        MPI_Allreduce(&(n_interactions), &globalChecksum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ThisTask==0) {
//             std::cout << "interaction checksum:" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << globalChecksum << " vs 1372670968 (64)" << std::endl;
            std::cout << "interaction checksum:" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << globalChecksum << " vs 260292 (4)" << std::endl;
        }
    }
    
    {
        long double localChecksum = 0.;
        for(i = 0; i < nActiveTot; i++) {
            localChecksum += local_AGN_alltau[i];
        }
        long double globalChecksum=0.;
        MPI_Allreduce(&(localChecksum), &globalChecksum, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (ThisTask==0) {
//             std::cout << "tau checksum:" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << globalChecksum << " vs 12592072.70693498738 (64)" << std::endl;
            std::cout << "tau checksum:" << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << globalChecksum << " vs 74.84413335002857746 (4)" << std::endl;
        }
    }
#endif

    tstart = system_time();
    MPI_Barrier(MPI_COMM_WORLD);
    tend = system_time();
    tdiffwait += dtime_system(tstart,tend);

    tstart = system_time();
    // for large numbers of processors, this might be slow
    // but currently it's way faster than the raytracing, so it doesn't really matter
    MPI_Allreduce(local_AGN_alltau, AGN_alltau, nActiveTot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tend = system_time();
    tdiffcomm += dtime_system(tstart,tend);

    

    // calculate column density to each active local particle
    // store this, and calculate heating and cooling and rad pressure from this in call from kicks.c
    jp = 0;
//     std::ofstream f;
//     f.open("depth_dump"+std::to_string(ThisTask)+".dat");

    for(i = 0; i < NumPart; i++) {
//        if( isGas(i) && isAlive(i) && TimeBinActive[P[i].TimeBin]  && isDusty(i)) {
#ifdef SOTON_AGN_DOUBLESTART
        if( isGas(i) && isAlive(i) && ( TimeBinActive[P[i].TimeBin] || All.NumCurrentTiStep<=1)) {
#else
        if( isGas(i) && isAlive(i) && isActive(i) ) {
#endif
            
            depths[i] = AGN_alltau[jp+localOffset];
//             f << P[i].ID << " " << AGN_alltau[jp+localOffset] << std::endl;
            jp++;
        }
    }
    
    delete[] local_AGN_alltau;
    delete[] AGN_alltau;
    delete[] AllPos;


    t1 = system_time();
    tdiffray = dtime_system(t0,t1)-tdiffbeam-tdiffcomm-tdiffwait;

    CPU_Step[CPU_SOTONAGNRAY] += tdiffray;
    CPU_Step[CPU_SOTONAGNBEAM] += tdiffbeam;
    CPU_Step[CPU_SOTONAGNCOMM] += tdiffcomm;
    CPU_Step[CPU_SOTONAGNWAIT] += tdiffwait;

     #ifdef BENCHMARK
         exit(0);
     #endif

}

