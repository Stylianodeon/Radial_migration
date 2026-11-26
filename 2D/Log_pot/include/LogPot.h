#pragma once
#include <string>
#include "VecUtils.h"

//global variables
extern double m;
extern double G;
extern double R_o;
extern double R_sigma;
extern double R_CR;
extern double Sigma_o;
extern double sigma_o;
extern double e_s;
extern double theta_deg;
extern double theta;
extern double alpha;
extern double omega_p;
extern double pert_strength;
extern double v_c;
extern double R_sp_taper;
extern double n_sp_taper;

// Function declarations
double total_potential_energy(const Vec2&, const double&);
double potential_energy_unperturbed(const Vec2&);
double total_PE_time(const Vec2&, const double&);
Vec2 LogPot_acc(const Vec2&);
void Leapfrog_integrator_unperturbed(Vec2&, Vec2&, const double&);
void Leapfrog_integrator_perturbed(Vec2&, Vec2&, const double&, const double&);
Vec2 total_acceleration(const Vec2&, const double&);

void simulate_trajectory(Vec2& pos, Vec2& vel, const std::string& basename, double& total_time, double& dt, bool use_pert);

void integrate_no_io(Vec2 pos, Vec2 vel, double total_time, double dt, bool use_pert, double& Lz0, double& Lz_end, double& JR0, double& JR_end,
                     double& eR0, double& eR_end, Vec2* pos_end = nullptr, Vec2* vel_end = nullptr);
