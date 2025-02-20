
/**
 * @brief Defines the following as free functions
 *	  - singular directions (Xi), 
 *	  - indicator functions (IN_BETAi) for each region
 *	  - integer functions used for bouncing trajectories (alpha, gamma, kappa)
 *	  - areas of each region (A_Bi)
 *
 * TODO: Refactor 
 */

#ifndef COORDINATESPACE_H_INCLUDED
#define COORDINATESPACE_H_INCLUDED

#include "config.h"
// ===================================================================================
// SINGULAR DIRECTIONS. START.
// ===================================================================================
constexpr float_ X0 = -PI2;
 
float_ X1(const float_& d, const float_& h) { return -atan(h/d); }
 
float_ X2(const float_& d, const float_& h) { return -atan((h+d)/(4*d)); }
 
float_ X3(const float_& d, const float_& h) { return -atan((h+d)/(5*d)); }
 
float_ X4(const float_& d, const float_& h) { return -atan((h+d-DELTA)/(d+DELTA)); }
 
float_ X5(const float_& d, const float_& h) { return -atan((h)/(d+1)); }
 
float_ X6(const float_& d, const float_& h) { return -atan((h)/(2*d+1)); }
 
float_ X7(const float_& d, const float_& h) { return atan((d-h)/(2*d+1)); }
 
float_ X8(const float_& d, const float_& h) { return atan((DELTA-h)/(2*d+DELTA)); }
 
float_ X9(const float_& d, const float_& h) { return atan((DELTA-h)/(d+DELTA)); }
 
float_ X10(const float_& d, const float_& h) { return atan((1-h)/(2*d+1)); }
 
float_ X11(const float_& d, const float_& h) { return atan((1-h)/(d+1)); }
 
float_ X12(const float_& d, const float_& h) { return atan((d-h)/d); }
 
float_ X13(const float_& d, const float_& h) { return atan((d+DELTA-h)/(d+DELTA)); }
 
float_ X14(const float_& d, const float_& h) { return atan((d+1-h)/(d+1)); }
 
float_ X15(const float_& d, const float_& h) { return atan(d+1-h); }
 
float_ X16(const float_& d, const float_& h) { return atan((d+DELTA-h)/DELTA); }
 
float_ X17(const float_& d, const float_& h) { return atan(2*d+1-h); }
 
float_ X18(const float_& d, const float_& h) { return atan((DELTA+2*d-h)/DELTA); }
constexpr float_ X19 = PI2;
// ===================================================================================
// SINGULAR DIRECTIONS. END.
// ===================================================================================
// ===================================================================================
// INDICATOR FUNCTIONS FOR EACH REGION. START.
// ===================================================================================
bool IN_BETA0(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X0 && theta < X1(d, h))); 
}
/// Individual bouncing trajectories not distinguished at this stage
bool IN_BETA1(const float_& d, const float_& h, const float_& theta) { 
    return ((h > (1/(float_)6) && h < d) && (theta > X1(d, h) && theta < X2(d, h)));
}
 
bool IN_Y1(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X1(d, h) && theta < X4(d, h))); 
}
 
bool IN_BETA2(const float_& d, const float_& h, const float_& theta) { 
    return (((h > 1/(float_)8 && h < d) && (theta > X2(d, h) && theta < X3(d, h))) && IN_Y1(d, h, theta)); 
}
 
bool IN_BETA3(const float_& d, const float_& h, const float_& theta) { 
    return (((h > 0 && h < (1/(float_)3)) && (theta > X3(d, h) && theta < X4(d, h))) && IN_Y1(d, h, theta)); 
}
 
bool IN_BETA4(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X4(d, h) && theta < X5(d, h))); 
}
 
bool IN_BETA5(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X5(d, h) && theta < X6(d, h))); 
}
 
bool IN_BETA6(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X6(d, h) && theta < X7(d, h))); 
}
 
bool IN_BETA7(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X7(d, h) && theta < X8(d, h))); 
}
 
bool IN_BETA8(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X8(d, h) && theta < X9(d, h))); 
}
 
bool IN_Y7(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X9(d, h) && theta < X12(d, h))); 
}
 
bool IN_BETA9(const float_& d, const float_& h, const float_& theta) { 
    return (((h > 0 && h < d) && (theta > X9(d, h) && theta < X10(d, h))) && IN_Y7(d,h,theta)); 
}
 
bool IN_BETA10(const float_& d, const float_& h, const float_& theta) { 
    return (((h > 0 && h < 1/(float_)3) && (theta > X10(d, h) && theta < X11(d, h))) && IN_Y7(d,h,theta)); 
}
 
bool IN_BETA11(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < 1/(float_)4) && (theta > X11(d, h) && theta < X12(d, h))); 
}
/// Individual bouncing ball trajectories not distinguished at this stage
bool IN_BETA12(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X12(d, h) && theta < X13(d, h))); 
}
 
bool IN_BETA13(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X13(d, h) && theta < X14(d, h))); 
}
 
bool IN_BETA14(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X14(d, h) && theta < X15(d, h))); 
}
 
bool IN_BETA15(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X15(d, h) && theta < X16(d, h))); 
}
 
bool IN_BETA16(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X16(d, h) && theta < X17(d, h))); 
}
 
bool IN_BETA17(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X17(d, h) && theta < X18(d, h))); 
}
 
bool IN_BETA18(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X18(d, h) && theta < X19)); 
}
/// The following regions are for when d = 1. These differ from above when first collision is on \Gamma_0
bool IN_A(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > -atan(h) && theta < atan((float_)1 - 2*h))); 
}
 
bool IN_ZJ0(const float_& d, const float_& h, const float_& theta) { // REGION 4
    return (((h > 0 && h < d) && (theta > X1(d, h) && theta < X5(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ1(const float_& d, const float_& h, const float_& theta) { // REGION 5
    return (((h > 0 && h < 2/(float_)3) && (theta > X5(d, h) && theta < X6(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ2(const float_& d, const float_& h, const float_& theta) { // REGION 6
    return (((h > 0 && h < 3/(float_)5) && (theta > X6(d, h) && theta < X8(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ3(const float_& d, const float_& h, const float_& theta) { // REGION 7
    return (((h > 0 && h < 1/(float_)2) && (theta > X8(d, h) && theta < X9(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ4(const float_& d, const float_& h, const float_& theta) { // REGION 8
    return (((h > 0 && h < 1/(float_)2) && (theta > X9(d, h) && theta < X7(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ5(const float_& d, const float_& h, const float_& theta) { // REGION 9
    return (((h > 0 && h < 2/(float_)5) && (theta > X7(d, h) && theta < X11(d,h))) && IN_A(d,h,theta)); 
}
 
bool IN_ZJ6(const float_& d, const float_& h, const float_& theta) { // REGION 10
    return (((h > 0 && h < 1/(float_)3) && (theta > X11(d, h) && theta < atan((DELTA-h)/DELTA))) && IN_A(d,h,theta)); 
}
 
bool IN_BETA19(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > atan((DELTA-h)/DELTA) && theta < atan(d-h))); 
}

/// Regions determining a first collision on \Gamma_0, \Gamma_2 or \Gamma_3
bool IN_G0(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X0 && theta < atan((DELTA-h)/DELTA))); 
}
 
bool IN_G2(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > atan((DELTA-h)/DELTA) && theta < X16(d,h))); 
} 
 
bool IN_G3(const float_& d, const float_& h, const float_& theta) { 
    return ((h > 0 && h < d) && (theta > X16(d,h) && theta < X19)); 
} 
// ===================================================================================
// INDICATOR FUNCTIONS FOR EACH REGION. END.
// ===================================================================================
// ===================================================================================
// INTEGER FUNCTIONS FOR BOUNCING TRAJECTORIES. START.
// ===================================================================================
/// Equation 6
int_ll alpha(const float_& d, const float_& h, const float_& theta) { 
    return ceil((-h - tan(theta)) / (d * (1 + tan(theta))));
}
/// Equation 11
int_ll kappa(const float_& d, const float_& h, const float_& theta) { 
    return ceil((d - h - tan(theta)) / (d*(tan(theta) - 1)));
}
/// Equation 13
int_ll gamma(const float_& d, const float_& h, const float_& theta) { 
    return ceil((h - 1 - d) / (d*(1 - tan(theta))));
}
// ===================================================================================
// INTEGER FUNCTIONS FOR BOUNCING TRAJECTORIES. END.
// ===================================================================================
// ===================================================================================
// REGION AREAS. START.
// ===================================================================================

/// \int atan(ax+b) dh = ...
template <typename T> 
T ataxb(const T& x, const T& a, const T& b) { 
    T axb = a*x + b;
    return ((2*(axb)*atan(axb) - std::log(axb*axb + 1))/(2*a));
}

 // \int_0^d (-X1+pi/2)dh
float_ A_B0(const float_& d) {  // ~0.56599 (d = 1/2) 1.1320 (d = 1)
    return d*(PI + std::log(4))/((float_)4); 
}

 // \int_0^d (X2 - X1) dh 
float_ A_B1(const float_& d) { // ~0.06076 (d = 1/2)
    float_ la = 1/(4*d);
    float_ lb = 1/(float_)4;
    float_ ra = 1/d;
    float_ rb = 0;
    return -ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	   - (-ataxb<float_>(1/(float_)6,la,lb) + ataxb<float_>(1/(float_)6,ra,rb)); 
}

 // \int_0^(1/3) (-X4 + X3) dh - \int_0^(1/8) (-X1 + X3) dh 
float_ A_B3(const float_& d) { // ~ 0.0200 (d = 1/2)
    float_ la = 1/(d+DELTA);
    float_ lb = (d-DELTA)/(d+DELTA);
    float_ ra = 1/(5*d);
    float_ rb = 1/(float_)5;
    float_ x4mx3 = -ataxb<float_>(1/(float_)3,la,lb) + ataxb<float_>(1/(float_)3,ra,rb) 
	           - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.03218
    la = 1/d;
    lb = 0;
    float_ x1mx3 = -ataxb<float_>(1/(float_)8,la,lb) + ataxb<float_>(1/(float_)8,ra,rb) 
	           - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.01219
    return x4mx3 - x1mx3; //0.019990
}

 // \int_0^d (-X4 + X1) dh - A_B3 - A_B1
float_ A_B2(const float_& d) { // ~0.01840 (d = 1/2)
    float_ la = 1/(d+DELTA);
    float_ lb = (d-DELTA)/(d+DELTA);
    float_ ra = 1/d;
    float_ rb = 0;
    float_ x4mx0 = -ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	           - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.09916
    return x4mx0 - A_B3(d) - A_B1(d); 
}

 // \int_0^d (-X5 + X4) dh 
float_ A_B4(const float_& d) { // ~0.03840 (d = 1/2)
    float_ la = 1/(d+1);
    float_ lb = 0;
    float_ ra = 1/(d+DELTA);
    float_ rb = (d-DELTA)/(d+DELTA);
    return -ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	   - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (-X6 + X5) dh 
float_ A_B5(const float_& d) { // ~0.0200 (d = 1/2)
    float_ la = 1/(2*d+1);
    float_ lb = 0;
    float_ ra = 1/(d+1);
    float_ rb = 0;
    return -ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	   - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X7 + X6) dh 
float_ A_B6(const float_& d) { // ~0.1237 (d = 1/2)
    float_ la = -1/(2*d+1);
    float_ lb = d/(2*d+1);
    float_ ra = 1/(2*d+1);
    float_ rb = 0;
    return ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X8 - X7) dh 
float_ A_B7(const float_& d) { // ~0.0200 (d = 1/2) 
    float_ la = -1/(2*d+DELTA);
    float_ lb = DELTA/(2*d+DELTA);
    float_ ra = -1/(2*d+1);
    float_ rb = d/(2*d+1);
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X9 - X8) dh 
float_ A_B8(const float_& d) { // ~0.03840 (d = 1/2) 
    float_ la = -1/(d+DELTA);
    float_ lb = DELTA/(d+DELTA);
    float_ ra = -1/(2*d+DELTA);
    float_ rb = DELTA/(2*d+DELTA);
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

// A_B9,10,11
 // \int_0^d (X10 - X9) dh - \int_(1/3)^d (X10-X12) dh
float_ A_B9(const float_& d) { // ~0.038397  (d = 1/2)
    float_ la = -1/(2*d+1);
    float_ lb = -la;
    float_ ra = -1/(d+DELTA);
    float_ rb = DELTA/(d+DELTA);
    float_ x10mx9 = ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	            - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); // 0.05839
    ra = -1/d;
    rb = 1;
    float_ x10mx12 = ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	             - (ataxb<float_>(1/(float_)3,la,lb) - ataxb<float_>(1/(float_)3,ra,rb)); //0.0200
    return x10mx9 - x10mx12; 
}

 // \int_0^(1/4) (X12 - X11) dh 
float_ A_B11(const float_& d) { // ~0.02746  (d = 1/2)
    float_ la = -1/d;
    float_ lb = 1;
    float_ ra = -1/(d+1);
    float_ rb = -ra;
    return ataxb<float_>(1/(float_)4,la,lb) - ataxb<float_>(1/(float_)4,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X12 - X9) dh - A_B11 - A_B9
float_ A_B10(const float_& d) { // ~0.0333  (d = 1/2)
    float_ la = -1/d;
    float_ lb = 1;
    float_ ra = -1/(d+DELTA);
    float_ rb = DELTA/(d+DELTA);
    float_ x12mx10 = ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	             - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
    return x12mx10 - A_B11(d) - A_B9(d);
}

 // \int_0^d (X13 - X12) dh 
float_ A_B12(const float_& d) { // ~0.09916 (d = 1/2), 0.13755 (d = 1)
    float_ la = -1/(d + DELTA);
    float_ ra = -1/d;
    return ataxb<float_>(d,la,1) - ataxb<float_>(d,ra,1) 
	   - (ataxb<float_>(0,la,1) - ataxb<float_>(0,ra,1)); 
}

 // \int_0^d (X14 - X13) dh 
float_ A_B13(const float_& d) { // ~0.02746 (d = 1/2), 0.06076 (d = 1)
    float_ la = -1/(d + 1);
    float_ ra = -1/(d + DELTA);
    return ataxb<float_>(d,la,1) - ataxb<float_>(d,ra,1) 
	   - (ataxb<float_>(0,la,1) - ataxb<float_>(0,ra,1)); 
}

 // \int_0^d (X15 - X14) dh 
float_ A_B14(const float_& d) { // ~ 0.10001 (d = 1/2), 0.3336 (d = 1)
    float_ lb = d+1;
    float_ ra = -1/(d+1);
    return ataxb<float_>(d,-1,lb) - ataxb<float_>(d,ra,1) 
	   - (ataxb<float_>(0,-1,lb) - ataxb<float_>(0,ra,1)); 
}

 // \int_0^d (X16 - X15) dh 
float_ A_B15(const float_& d) { // ~0.039334 (d = 1/2), 0.10775 (d = 1)
    float_ la = -1/DELTA;
    float_ lb = (d+DELTA)/DELTA;
    float_ ra = -1;
    float_ rb = d+1;
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X17 - X16) dh 
float_ A_B16(const float_& d) { // ~0.039338 (d = 1/2), 0.10775 (d = 1)
    float_ la = -1;
    float_ lb = 2*d + 1;
    float_ ra = -1/DELTA;
    float_ rb = (d+DELTA)/DELTA;
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X18 - X17) dh 
float_ A_B17(const float_& d) { // ~0.06842 (d = 1/2), 0.13479 (d = 1)
    float_ la = -1/DELTA;
    float_ lb = (DELTA+2*d)/DELTA; 
    float_ ra = -1;
    float_ rb = 2*d + 1;
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (pi/2 - X18) dh 
float_ A_B18(const float_& d) { // ~0.19227 (d = 1/2), 0.2497 (d = 1)
    float_ ra = -1/DELTA;
    float_ rb = (DELTA+2*d)/DELTA; 
    return d*PI2 - ataxb<float_>(d,ra,rb) + ataxb<float_>(0,ra,rb); 
}


 // int_0^d (-X4 + X0) dh - \int_(2/3)^d (-X4 - X19) dh
float_ A_ZJ0(const float_& d) { // ~0.13756 (d = 1)
    float_ la = 1/(d+1);
    float_ lb = 0;
    float_ ra = 1/d;
    float_ rb = 0;
    float_ x4x0 = -ataxb<float_>(d,la,lb) + ataxb<float_>(d,ra,rb) 
	          - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.19832
    ra = -1/DELTA;
    rb = 1;
    float_ x4x19 = -ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	           - (-ataxb<float_>(2/(float_)3,la,lb) - ataxb<float_>(2/(float_)3,ra,rb)); //~0.06076
    return x4x0 - x4x19;
}

 // int_0^(3/5) (-X6 + X5) dh + \int_(3/5)^(2/3) (X19 + X5) dh
float_ A_ZJ1(const float_& d) { // ~0.03219 (d = 1)
    float_ la = 1/(2*d+1);
    float_ lb = 0;
    float_ ra = 1/(d+1);
    float_ rb = 0;
    float_ x6x5 = -ataxb<float_>(3/(float_)5,la,lb) + ataxb<float_>(3/(float_)5,ra,rb) 
	          - (-ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.04113
    la = -1/DELTA;
    lb = 1;
    float_ x19x5 = ataxb<float_>(2/(float_)3,la,lb) + ataxb<float_>(2/(float_)3,ra,rb) 
		   - (ataxb<float_>(3/(float_)5,la,lb) + ataxb<float_>(3/(float_)5,ra,rb)); //~0.0031
    return x6x5 + x19x5;
}

 // int_0^(1/2) (X8 + X6) dh + \int_(1/2)^(3/5) (X19 + X6) dh
float_ A_ZJ2(const float_& d) { // ~0.04422 (d = 1)
    float_ la = -1/(2*d + DELTA);
    float_ lb = DELTA/(2*d + DELTA);
    float_ ra = 1/(2*d+1);
    float_ rb = 0;
    float_ x8x6 = ataxb<float_>(1/(float_)2,la,lb) + ataxb<float_>(1/(float_)2,ra,rb) 
		  - (ataxb<float_>(0,la,lb) + ataxb<float_>(0,ra,rb)); //~0.09115
    la = -1/DELTA;
    lb = 1;
    float_ x19x6 = ataxb<float_>(3/(float_)5,la,lb) + ataxb<float_>(3/(float_)5,ra,rb) 
	           - (ataxb<float_>(1/(float_)2,la,lb) + ataxb<float_>(1/(float_)2,ra,rb)); //~0.00820
    return x8x6 + x19x6;
}

 // \int_0^(1/2) (X9 - X8) dh 
float_ A_ZJ3(const float_& d) { // ~0.032183 (d = 1)
    float_ la = -1/(d+DELTA);
    float_ lb = DELTA/(d+DELTA);
    float_ ra = -1/(2*d + DELTA);
    float_ rb = DELTA/(2*d + DELTA);
    return ataxb<float_>(1/(float_)2,la,lb) - ataxb<float_>(1/(float_)2,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb));
}

 // int_0^(2/5) (X7 - X9) dh + \int_(2/5)^(1/2) (X19 - X9) dh
float_ A_ZJ4(const float_& d) { // ~0.032183 (d = 1) 
    float_ la = -1/(2*d+1);
    float_ lb = d/(2*d+1);
    float_ ra = -1/(d + DELTA);
    float_ rb = DELTA/(d + DELTA);
    float_ x7x9 = ataxb<float_>(2/(float_)5,la,lb) - ataxb<float_>(2/(float_)5,ra,rb) 
	          - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
    la = -1/DELTA;
    lb = 1;
    float_ x19x9 = ataxb<float_>(1/(float_)2,la,lb) - ataxb<float_>(1/(float_)2,ra,rb) 
	           - (ataxb<float_>(2/(float_)5,la,lb) - ataxb<float_>(2/(float_)5,ra,rb)); //
    return x7x9 + x19x9;
}

 // int_0^(1/3) (X11 - X7) dh + \int_(1/3)^(2/5) (X19 - X7) dh
float_ A_ZJ5(const float_& d) { // ~0.04461 (d = 1)
    float_ la = -1/(d+1);
    float_ lb = 1/(d+1);
    float_ ra = -1/(1+2*d);
    float_ rb = d/(1+2*d);
    float_ x11x7 = ataxb<float_>(1/(float_)3,la,lb) - ataxb<float_>(1/(float_)3,ra,rb) 
		   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); //~0.04113
    la = -1/DELTA;
    lb = 1;
    float_ x19x7 = ataxb<float_>(2/(float_)5,la,lb) - ataxb<float_>(2/(float_)5,ra,rb) 
	           - (ataxb<float_>(1/(float_)3,la,lb) - ataxb<float_>(1/(float_)3,ra,rb)); //~0.00348
    return x11x7 + x19x7;
}

 // \int_0^(1/3) (X - X11) dh 
float_ A_ZJ6(const float_& d) { // ~0.0608  (d = 1)
    float_ la = -1/DELTA;
    float_ lb = 1;
    float_ ra = -1/(d+1);
    float_ rb = 1/(d+1);
    return ataxb<float_>(1/(float_)3,la,lb) - ataxb<float_>(1/(float_)3,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb)); 
}

 // \int_0^d (X19 - X18) dh 
float_ A_B19(const float_& d) { // ~0.43882 (d = 1)
    float_ la = -1;
    float_ lb = d;
    float_ ra = -1/DELTA;
    float_ rb = 1;
    return ataxb<float_>(d,la,lb) - ataxb<float_>(d,ra,rb) 
	   - (ataxb<float_>(0,la,lb) - ataxb<float_>(0,ra,rb));
}
// ===================================================================================
// REGION AREAS. END.
// ===================================================================================

#endif
