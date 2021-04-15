/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   intrinsic_depth.h
 * Author: davidjwilliamson
 *
 * Created on 14 April 2021, 11:31
 */

#ifndef INTRINSIC_DEPTH_H
#define INTRINSIC_DEPTH_H

#include "raytrace_helpers.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef INTRINSIC_DEPTH
double intrinsic_depth_angle(double theta) {
    if ( theta<All.intrinsic_tau_transition ) {
        return All.intrinsic_tau0-theta/All.intrinsic_index0;
    } else {
        return All.intrinsic_tau0-All.intrinsic_tau_transition/All.intrinsic_index0-(theta-All.intrinsic_tau_transition)/All.intrinsic_index1;
    }
}


double intrinsic_depth(double *pos) {
    double r2d = sqrt(square(pos[0])+square(pos[1]));
    double z = fabs(pos[2]);
    double theta = atan2(z,r2d);
    return intrinsic_depth_angle(theta);
}
#else

double inline intrinsic_depth_angle(double theta) {
    return 0;
}

double inline intrinsic_depth(double *pos) {
    return 0;
}
#endif



#ifdef __cplusplus
}
#endif

#endif /* INTRINSIC_DEPTH_H */

