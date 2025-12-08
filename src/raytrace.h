#pragma once
#include <glm/glm.hpp>

// Dneg wormhole without gravity, metric:
// ds^2 = -dt^2 + dℓ^2 + r(ℓ)^2 (dθ^2 + sin^2θ dφ^2)
//
// r(ℓ) is defined piecewise: cylindrical interior, smoothed transition outside.

struct WormholeParams {
    float rho;   // throat radius
    float a;     // half-length of cylindrical interior
    float M;     // transition "mass" parameter (controls lensing width)
};

// A single ray's state:
struct RayState {
    float l;       // ℓ
    float theta;   // θ
    float phi;     // φ
    float p_l;     // p_ℓ
    float p_theta; // p_θ
    float b;       // p_φ (impact parameter around polar axis)
};

// A struct that stores ray positions
struct RayPosition {
    glm::vec4 pos; // the Euclidean coords in world space
    float l;
};

// Compute r(l) as in Equation 5a, 5b, 5c
float r_of_l(float l, const WormholeParams &p);

// Compute dr/dl as in Equation 5a
float drdl_of_l(float l, const WormholeParams &p);

// Compute derivatives d/dt of the state as in Equation A.7
void geodesicDeriv(const RayState &s,
                   const WormholeParams &wp,
                   float &d_l,
                   float &d_theta,
                   float &d_phi,
                   float &d_p_l,
                   float &d_p_theta);

// One Euler's method step for numerically solving the geodesic ODE
void odeStep(RayState &s, const WormholeParams &wp, float dt);

// One RK4 step for numerically solving the geodesic ODE
void rk4Step(RayState &s, const WormholeParams &wp, float dt);
