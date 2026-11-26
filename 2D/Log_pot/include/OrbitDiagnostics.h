#pragma once
#include "VecUtils.h"

double angular_momentum(const Vec2& pos, const Vec2& vel);
double kinetic_energy(const Vec2& vel);
double jacobi_integral(const Vec2& pos, const Vec2& vel, const double& sim_time, const double& total_time);
double getJRunperturbed(Vec2 pos, Vec2 vel, double dt);
double random_energy(const Vec2& pos, const Vec2& vel);
