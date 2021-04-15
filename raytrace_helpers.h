#ifndef RAYTRACE_HELPERS_H
#define RAYTRACE_HELPERS_H


// helper functions
static inline double square(double x) {
    return x*x;
}

static inline double cube(double x) {
    return x*x*x;
}


static inline double powerfour(double x) {
    return square(square(x));
}


static inline int square(int x) {
    return x*x;
}

static inline double get_d12_norm2(double *r1, double *r2) {
    double norm=0.;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        norm+=square(r1[ii]-r2[ii]);
    }
    
    return norm;
}

static inline double get_norm2(double *r1) {
    double norm=0.;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        norm+=square(r1[ii]);
    }
    
    return norm;
}

// return values:
// 0 = not between
// 1 = strictly between (entire kernel is between centres)
// 2 = partially between (kernel overlaps one side)
static inline int between_particles(double *r1,double *r2,double *r3,double h_3, double r12) {


// check pointers do the correct thing
    double z1=0.,z2=0.,z3=0.;

    for ( int jj=0 ; jj<3 ; jj++ ) {
        double dd = r1[jj]-r2[jj];
        z1+=r1[jj]*dd;
        z2+=r2[jj]*dd;
        z3+=r3[jj]*dd;
    }
    
    z1/=r12;
    z2/=r12;
    z3/=r12;
    
    if ( z1>z2 ) { // source is more positive than target

        if ( z1>z3+h_3 && z3-h_3>z2 ) {
            return 1;
        } else if ( z1>z3-h_3 && z3+h_3>z2 ) {
            return 2;
        }
    } else if ( z1<z2 ) { // target is more positive from origin than source

        if ( z1<z3-h_3 && z3+h_3<z2 ) {
            return 1;
        } else if ( z1<z3+h_3 && z3-h_3<z2 ) {
            return 2;
        }
    } else {
        // particles are at a tangent from the source
        // part of the particle really should be providing extinction
        // but we're ignoring that for now
        return 0;
    }
    return 0;
}
 
// return values:
// 0 = not between
// 1 = strictly between (entire kernel is between centres)
// 2 = partially between (kernel overlaps one side)
// - same as between_particles, but assume r1=0
static inline int between_particle_agn(double *r2,double *r3,double h_3, double r) {
    double z3=0.; //z1=0, constant

    for ( int jj=0 ; jj<3 ; jj++ ) {
        z3+=-r3[jj]*r2[jj];
    }
    
    double z2 = -r;
    
    z3/=r;
    
    if ( z2<0 ) { // source is more positive than target

        if ( z3+h_3<0 && z3-h_3>z2 ) {
            return 1;
        } else if ( z3-h_3<0 && z3+h_3>z2 ) {
            return 2;
        }
    } else if ( z2>0 ) { // target is more positive from origin than source

        if ( z3-h_3>0 && z3+h_3<z2 ) {
            return 1;
        } else if ( z3+h_3>0 && z3-h_3<z2 ) {
            return 2;
        }
    } else {
        // particles are at a tangent from the source
        // part of the particle really should be providing extinction
        // but we're ignoring that for now
        return 0;
    }
    return 0;
}
 

static inline double intersect_agn_d2_nonorm(double *r2,double *r3) {
    double dr32[3];
    double crossp;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        dr32[ii] = r3[ii]-r2[ii];
    }
    
    crossp = square(   r3[1]*dr32[2]-  r3[2]*dr32[1]);
    crossp+= square(-  r3[0]*dr32[2]+  r3[2]*dr32[0]);
    crossp+= square(   r3[0]*dr32[1]-  r3[1]*dr32[0]);
    
    return crossp;
}


static inline double intersect_d2_nonorm(double *r1,double *r2,double *r3) {
    double dr31[3],dr32[3];
    double crossp;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        dr31[ii] = r3[ii]-r1[ii];
        dr32[ii] = r3[ii]-r2[ii];
    }
    
    crossp = square( dr31[1]*dr32[2]-dr31[2]*dr32[1]);
    crossp+= square(-dr31[0]*dr32[2]+dr31[2]*dr32[0]);
    crossp+= square( dr31[0]*dr32[1]-dr31[1]*dr32[0]);
    
    return crossp;
}
#endif