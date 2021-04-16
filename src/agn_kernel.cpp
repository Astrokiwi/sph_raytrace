#include "rt_prototypes.h"
#include <vector>
//#include <stdlib.h>

#include "coupling.h"
#include <math.h>
#include <iostream>

//#include "../../allvars.h"
//#include "../../proto.h"
//#include "../../kernel.h"

//#include <fstream>

// helper functions
static inline double square(double x) {
    return x*x;
}

static inline int square(int x) {
    return x*x;
}

void AGN_Kernel::integrate_flat() {

    flattened_table2.resize(L_flat);

    fat_table2.resize(L_flat, std::vector<double>(zL));

    
    double r3a,r3b;
    double r2;
    double rsq;
    
    double dz = 1./zL;
    
    for ( int ir = 0 ; ir<L_flat ; ir++ ) {
        flattened_table2[ir] = 0.;
        r2 = ((double)(ir))/L_flat;
        rsq = square(r2);
        for ( int iz = 0 ; iz<zL ; iz++ ) {
            //r3 = sqrt(square(((double)iz)/zL)+square(r2));
            //flattened_table[ir]+=2.*w(r3)*dz; // even around 0.

            r3a = sqrt(square(((double)iz)/zL)+rsq);
            r3b = sqrt(square(((double)(iz+1))/zL)+rsq);
            flattened_table2[ir]+=2.*(w(r3a)+w(r3b))/2.*dz; // even around 0.
            fat_table2[ir][iz] = flattened_table2[ir]/2.; // half sized table
        }
    }
    
}

void AGN_Kernel::set_L_flat(int L_in, int zL_in) {
    L_flat = L_in;
    zL = zL_in;
    integrate_flat();
}

AGN_Kernel::AGN_Kernel() {
    set_L_flat(256,4096); // default
}


        
double AGN_Kernel::w(double r, double h) {
    return kernel_wk(r,h);
}


        
double AGN_Kernel::w(double x) {
    return this->w(x,1.);
}

// skin_length is h-skin_depth
double AGN_Kernel::skin_rad(double h_in, double opac) {
    int ih;
    double fh;
    double logh;
    double h = h_in/opac; // effective h, see calc_depth_optical_table(...) below
    
    logh = log10(h);
    fh = this->nh*(logh-this->loghmin)/(this->loghmax-this->loghmin);
    
    ih = floor(fh);
    
    if ( ih>=this->nh-1 ) {
        return -1.;
    } else {
        fh = fh-ih;
        //std::cout << ih << " " << fh << " " << skin_table[ih] << " " << skin_table[ih+1] << std::endl;
        
        return (skin_table[ih]*(1.-fh)+skin_table[ih+1]*fh)*h;
    }
}

void AGN_Kernel::calc_depth_optical_table(double depth, double mass, double loghmin, double loghmax, int nh) {
/* This used tabulated by h for a fixed opacity.
But we can use this for variable opacity, because result depends only on h^2/mass/kappa (where kappa is opacity in area/mass units).
So we set kappa=1 in internal units for this run ( (1 kpc)**2/(1e10 Msun) = 0.48 cm**2/g ), and then the "effective h" that we use
to look up the skindepth table is just: h/sqrt(kappa).

This means we need a larger range of possible h values than before, but that's no problem.
*/
    double depthtarget = depth/mass/1.;
    
    int iz,jz;
    
    double low_tab,high_tab,fz;
    double realz;
    double h,h2,w1,wtarget;
    
    int ih;
    
    //dh = (loghmax-loghmin)/nh;
    
    this->loghmin = loghmin;
    this->loghmax = loghmax;
    this->nh = nh;
    skin_table.resize(nh);
    
    //for ( double logh=-6.; logh<-1. ; logh+=.1 ) {
    for ( ih=0 ; ih<nh ; ih++ ) {
    
        h=pow(10.,loghmin+(ih*(loghmax-loghmin))/nh);
    
        h2 = square(h);
        w1 = flattened_table2[0]/2.;
        wtarget = w1-depthtarget*h2;
        
        if ( wtarget<0. ) {
            realz = 0.;
        } else {
            
            // linear search - change to binary later
            jz=-1;
            for ( iz=0 ; iz<zL ; iz++ ) {
                if ( fat_table2[0][iz]>wtarget ) {
                    jz = iz-1;
                    iz = zL;
                }
            }
        
            if ( jz==-1 ) {
                low_tab = 0.;
            } else {
                low_tab = fat_table2[0][jz];
            }
        
            high_tab = fat_table2[0][jz+1];
        
            fz = (wtarget-low_tab)/(high_tab-low_tab);
        
    
            realz = ((jz+1+fz)*h)/zL;
        }
        
        skin_table[ih] = realz/h;
    }

}


double AGN_Kernel::flat_w2(double r2, double h2) {
    // from table
    double x = L_flat*r2/h2;
    int ix = (int)(x);
    if ( ix>=L_flat ) {
        return 0.;
    }
    
    // interpolate
    double klow,khigh,kf;
    klow = flattened_table2[ix];
    kf = x-ix;
    if ( ix==L_flat-1 ) {
        khigh = 0.;
    } else {
        khigh = flattened_table2[ix+1];
    }
    
    double k = klow*(1.-kf)+khigh*kf;
    if ( kf<0. || kf>1. || k<0. ) {
        std::cout << "wtf " << klow << " " << khigh << " " << kf << " " << x << " " << ix << std::endl;
        exit(1);
    }
    
    k/=h2; // flattened kernel has h=1, need to take real h into account
    return k;
}



double AGN_Kernel::half_flat_w2_truncated(double r2, double h2, double z) {
    if ( z==0. ) return 0.;
    if ( z<0. ) {
        return -this->half_flat_w2_truncated(r2,h2,-z);
    }
    // from table
    double x = L_flat*r2/h2;
    int ix = (int)(x);
    if ( ix>=L_flat ) {
        return 0.;
    }

    double zt = z*zL; 
    int iz = (int)(zt)-1;// trapezium integration - iz=0 is dz from centre plane
    if ( iz>=zL-1 ) {
        return this->flat_w2(r2,h2)/2.; // for hemisphere
    }
    
    double xf = x-ix;


    double klow,khigh;
    
    // interpolate - low z
    if ( iz<0 ) {
        klow = 0.;
        khigh = 0.;
    } else {
        klow = fat_table2[ix][iz];
        if ( ix==L_flat-1 ) {
            khigh = 0.;
        } else {
            khigh = fat_table2[ix+1][iz];
        }
    }
    
    double k0 = klow*(1.-xf)+khigh*xf;
    if ( xf<0. || xf>1. || k0<0. ) {
        std::cout << "wtf " << klow << " " << khigh << " " << xf << " " << x << " " << ix << std::endl;
        exit(1);
    }

    
    // interpolate - high z
    //double klow,khigh,kf;
    klow = fat_table2[ix][iz+1];
    if ( ix==L_flat-1 ) {
        khigh = 0.;
    } else {
        khigh = fat_table2[ix+1][iz+1];
    }
    
    double k1 = klow*(1.-xf)+khigh*xf;
    if ( xf<0. || xf>1. || k1<0. ) {
        std::cout << "wtf " << klow << " " << khigh << " " << xf << " " << x << " " << ix << std::endl;
        exit(1);
    }


    double zf = zt-(iz+1);
    double k = k0*(1.-zf)+k1*zf;
    
    k/=h2; // flattened kernel has h=1, need to take real h into account
    return k;
}
