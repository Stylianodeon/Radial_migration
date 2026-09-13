#pragma once
#include <vector>
#include <tuple>
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
double effective_potential(const Vec2&, const bool&);
Vec2 effective_potential_gradient(const Vec2&, const bool&);
bool is_near_cluster(const std::vector<std::tuple<double, double, double>>&, const std::tuple<double, double, double>&, double);
void scan_Lagrange_points(double, double, double, double, double, const std::string&, bool);
double perturbation_scaling_factor(double);
double spiral_radial_envelope(double);
double d_spiral_radial_envelope_dR(double);
double total_potential_energy(const Vec2&, const double&);
double potential_energy_unperturbed(const Vec2&);
double total_PE_time(const Vec2&, const double&);
Vec2 LogPot_acc(const Vec2&);
void Leapfrog_integrator_unperturbed(Vec2&, Vec2&, const double&);
void Leapfrog_integrator_perturbed(Vec2&, Vec2&, const double&, const double&);
double angular_momentum(const Vec2&, const Vec2&);
double kinetic_energy(const Vec2&);
Vec2 total_acceleration(const Vec2&, const double&);
double jacobi_integral(const Vec2&, const Vec2&, const double&, const double&);
double effective_potential_time(const Vec2&, const bool&, double);
Vec2   effective_potential_gradient_time(const Vec2&, const bool&, double, double);


int run_logpot_tools(int argc, char** argv);



/*

//Adding a perturbation with constant pattern speed
// Calculate accelerations using formula acceleration = - ∇Φ0 - ∇Φ1 - coriolis - centtrifugal force 
Vec3 coriolis_acceleration(const Vec3& cart_vel)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p};
    Vec3 coriolis_acc = cross_product(omega_p_cart, cart_vel);

    return coriolis_acc * -2.0;
}

Vec3 centrifugal_acceleration(const Vec3 &cart_pos)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p};
    Vec3 intermediate = cross_product(omega_p_cart, cart_pos);
    Vec3 centrifugal_acc = cross_product(omega_p_cart, intermediate);
    return centrifugal_acc * -1.0 ;
}


Vec3 perturbation_acceleration(const Vec3 &cart_pos)
{   
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);
    //Determine the potential amplitude and its derivative from Daniel & Wyse (2015) eqn. 37

    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R/R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
    double dPhi_S_dR = (2.0 * M_PI * G * e_s * Sigma * (cyl_pos.R-R_o)) / (R_o * m * (1.0/tan(theta)));

    //Determine partial derivatives of potential perturbation Φ1 in the rotating fram give from Dnaiel & Wyse (2015) eqn. 6
    double dPhi_dR = dPhi_S_dR * cos(m * cyl_pos.phi);
    double dPhi_dphi = -m * Phi_S * sin(m * cyl_pos.phi);
    //Determine derivatives for chain rule
    double dR_dx = cart_pos.x / cyl_pos.R;
    double dR_dy = cart_pos.y / cyl_pos.R;
    double dphi_dx = -cart_pos.y/(cyl_pos.R * cyl_pos.R);
    double dphi_dy = cart_pos.x /(cyl_pos.R * cyl_pos.R);
    //Calculate potential in cartesian coordinates
    double dPhi_dx = (dPhi_dR * dR_dx) + (dPhi_dphi * dphi_dx);
    double dPhi_dy = (dPhi_dR * dR_dy) + (dPhi_dphi * dphi_dy);

    return {0.0, 0.0, 0.0};      //{ -dPhi_dx, -dPhi_dy, 0.0 };
}

*/