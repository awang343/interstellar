#include "raytrace.h"
#include <cmath>

// Compute r(l) as in Equation 5a, 5b, 5c
float r_of_l(float l, const WormholeParams &p) {
    float absL = std::fabs(l);
    // Equation 5a
    if (absL <= p.a) {
        return p.rho;
    }

    // set up x = 2(|l| - a) / (π M)
    float x = 2.0 * (absL - p.a) / (M_PI * p.M);

    // r = ρ + M [ x arctan(x) - 0.5 ln(1 + x^2) ]
    float term = x * atan(x) - 0.5 * log(1.0 + x * x);
    return p.rho + p.M * term;
}

// Compute dr/dl as in Equation 5a
float drdl_of_l(float l, const WormholeParams &p) {
    float absL = std::fabs(l);
    if (absL <= p.a) {
        return 0.0;
    }

    // derivative wrt |l|-a of the integral form in 5a gives
    // d/d(|l|-a) r = (2/π) arctan(2(|l|-a)/(πM))
    float x = 2.0 * (absL - p.a) / (M_PI * p.M);
    float dterm_dabs = (2.0 / M_PI) * atan(x);

    float sign = (l >= 0.0) ? 1.0 : -1.0;
    return dterm_dabs * sign;
}

// Compute derivatives d/dt of the state as in Equation A.7
void geodesicDeriv(const RayState &s,
                   const WormholeParams &wp,
                   float &d_l,
                   float &d_theta,
                   float &d_phi,
                   float &d_p_l,
                   float &d_p_theta) {
    float r = r_of_l(s.l, wp);
    float drdl = drdl_of_l(s.l, wp);

    float sinT = sin(s.theta);
    float cosT = cos(s.theta);
    // avoid dividing by zero
    const float eps = 1e-6;
    if (std::fabs(sinT) < eps) {
        sinT = (sinT >= 0.0 ? eps : -eps);
    }

    // Compute B^2 = p_θ^2 + b^2 / sin^2θ as in Equation A.5b
    float B2 = s.p_theta * s.p_theta + (s.b * s.b) / (sinT * sinT);

    // Equations of motion
    d_l      = s.p_l;
    d_theta  = s.p_theta / (r * r);
    d_phi    = s.b / (r * r * sinT * sinT);
    d_p_l    = B2 * drdl / (r * r * r);
    d_p_theta= (s.b * s.b / (r * r)) * (cosT / (sinT * sinT * sinT));
}

// One Euler's method step for numerically solving the geodesic ODE
void odeStep(RayState &s, const WormholeParams &wp, float dt) {
    float k_l, k_theta, k_phi, k_p_l, k_p_theta;

    // obtain the partial derivatives
    geodesicDeriv(s, wp, k_l, k_theta, k_phi, k_p_l, k_p_theta);

    // update the state
    s.l += dt * k_l;
    s.theta += dt * k_theta;
    s.phi += dt * k_phi;
    s.p_l += dt * k_p_l;
    s.p_theta += dt * k_p_theta;
}

// One RK4 step for numerically solving the geodesic ODE
void rk4Step(RayState &s, const WormholeParams &wp, float dt) {
    float k1_l, k1_theta, k1_phi, k1_p_l, k1_p_theta;
    float k2_l, k2_theta, k2_phi, k2_p_l, k2_p_theta;
    float k3_l, k3_theta, k3_phi, k3_p_l, k3_p_theta;
    float k4_l, k4_theta, k4_phi, k4_p_l, k4_p_theta;

    // k1
    geodesicDeriv(s, wp, k1_l, k1_theta, k1_phi, k1_p_l, k1_p_theta);

    // k2
    {
        RayState tmp = s;
        tmp.l      += 0.5 * dt * k1_l;
        tmp.theta  += 0.5 * dt * k1_theta;
        tmp.phi    += 0.5 * dt * k1_phi;
        tmp.p_l    += 0.5 * dt * k1_p_l;
        tmp.p_theta+= 0.5 * dt * k1_p_theta;
        geodesicDeriv(tmp, wp, k2_l, k2_theta, k2_phi, k2_p_l, k2_p_theta);
    }

    // k3
    {
        RayState tmp = s;
        tmp.l      += 0.5 * dt * k2_l;
        tmp.theta  += 0.5 * dt * k2_theta;
        tmp.phi    += 0.5 * dt * k2_phi;
        tmp.p_l    += 0.5 * dt * k2_p_l;
        tmp.p_theta+= 0.5 * dt * k2_p_theta;
        geodesicDeriv(tmp, wp, k3_l, k3_theta, k3_phi, k3_p_l, k3_p_theta);
    }

    // k4
    {
        RayState tmp = s;
        tmp.l      += dt * k3_l;
        tmp.theta  += dt * k3_theta;
        tmp.phi    += dt * k3_phi;
        tmp.p_l    += dt * k3_p_l;
        tmp.p_theta+= dt * k3_p_theta;
        geodesicDeriv(tmp, wp, k4_l, k4_theta, k4_phi, k4_p_l, k4_p_theta);
    }

    s.l      += dt * (k1_l      + 2.0 * k2_l      + 2.0 * k3_l      + k4_l)      / 6.0;
    s.theta  += dt * (k1_theta  + 2.0 * k2_theta  + 2.0 * k3_theta  + k4_theta)  / 6.0;
    s.phi    += dt * (k1_phi    + 2.0 * k2_phi    + 2.0 * k3_phi    + k4_phi)    / 6.0;
    s.p_l    += dt * (k1_p_l    + 2.0 * k2_p_l    + 2.0 * k3_p_l    + k4_p_l)    / 6.0;
    s.p_theta+= dt * (k1_p_theta+ 2.0 * k2_p_theta+ 2.0 * k3_p_theta+ k4_p_theta)/ 6.0;
}
