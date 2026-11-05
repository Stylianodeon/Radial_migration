#pragma once

#include "VecUtils.h"
#include <string>

// Full-strength pattern at t=0, including the centrifugal potential.
double effective_potential(const Vec2& pos, const bool& use_pert);
Vec2 effective_potential_gradient(const Vec2& pos, const bool& use_pert);

// Time-dependent pattern, including the existing perturbation ramp.
double effective_potential_time(const Vec2& pos, const bool& use_pert, double sim_time);
Vec2 effective_potential_gradient_time(const Vec2& pos, const bool& use_pert, double sim_time, double total_time);

// Plot the gradient and solve instantaneous equilibria at the supplied snapshot time.
// Coordinates are inertial x,y. A zero-amplitude pattern is reported as a ring
// in comments, with no isolated point rows.
// Assumes finite, ordered bounds, a positive step, and valid model parameters
// (including a positive integer arm count m). Throws on numerical or output failure.
void scan_Lagrange_points(double xmin, double xmax, double ymin, double ymax, double step, const std::string& filename, bool use_pert, double sim_time);
