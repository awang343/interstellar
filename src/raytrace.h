#pragma once

// Dneg wormhole without gravity, metric:
// ds^2 = -dt^2 + dℓ^2 + r(ℓ)^2 (dθ^2 + sin^2θ dφ^2)
//
// r(ℓ) is defined piecewise: cylindrical interior, smoothed transition outside.
struct WormholeParams {
    double rho;   // throat radius
    double a;     // half-length of cylindrical interior
    double M;     // transition "mass" parameter (controls lensing width)
};

// A single ray's state:
struct RayState {
    double l;       // ℓ
    double theta;   // θ
    double phi;     // φ
    double p_l;     // p_ℓ
    double p_theta; // p_θ
    double b;       // p_φ (impact parameter around polar axis)
};

// Compute r(l) as in Equation 5a, 5b, 5c
double r_of_l(double l, const WormholeParams &p);

// Compute dr/dl as in Equation 5a
double drdl_of_l(double l, const WormholeParams &p);

// Compute derivatives d/dt of the state as in Equation A.7
void geodesicDeriv(const RayState &s,
                   const WormholeParams &wp,
                   double &d_l,
                   double &d_theta,
                   double &d_phi,
                   double &d_p_l,
                   double &d_p_theta);

// One Euler's method step for numerically solving the geodesic ODE
void odeStep(RayState &s, const WormholeParams &wp, double dt);

// One RK4 step for numerically solving the geodesic ODE
void rk4Step(RayState &s, const WormholeParams &wp, double dt);
