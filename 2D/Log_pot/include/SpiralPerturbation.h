#pragma once
#include "VecUtils.h"

double spiral_radial_envelope(double R);
double d_spiral_radial_envelope_dR(double R);
double perturbation_scaling_factor(double sim_time);

// A includes pert_strength and the supplied scale; dA is its radial derivative.
void spiral_amplitude(double R, double scale, double& A, double& dA);
double spiral_phase(double R, double phi, double time);
// Pass scale=1 for full amplitude, or perturbation_scaling_factor(time) for the ramp.
double spiral_potential(const Cyl& pos, double time, double scale);
Vec2 spiral_gradient(const Vec2& pos, double time, double scale);
