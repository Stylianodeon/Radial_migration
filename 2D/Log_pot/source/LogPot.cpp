#include <algorithm>
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <functional>
#include "VecUtils.h"
#include "LogPot.h"
#include "SpiralPerturbation.h"

using std::string;
using std::vector;
using std::ofstream;
using std::cos;
using std::sin;
using std::log;

// Global parameters
double m           = 4.0;
double G           = 4.5e-12;       // kpc^3 / (Msun * Myr^2)
double R_o         = 2.5;           // kpc
double R_sigma     = 3 * R_o;       // kpc
double R_CR        = 10;            // kpc
double Sigma_o     = 1.226e+9;      // Msun / kpc^2
double sigma_o     = 0.0;           // σ_R / 2 from Chiba et al. 2021
double e_s         = 0.1;           // old working value 0.3
double theta_deg   = 15;            // old working value 35
double theta       = theta_deg * M_PI / 180.0;
double alpha       = m / std::tan(theta);
double omega_p     = 0.022;         // pattern speed (rad/Myr)
double pert_strength = 0.1;
double v_c         = 0.22;          // kpc/Myr = 220 km/s

// Smooth radial envelope for the spiral perturbation only.
// W(R_sp_taper)=0.5. Larger n_sp_taper gives a sharper, but still smooth, cutoff.
double R_sp_taper = 5.0;            // kpc; tune this to suppress the low-Lz/inner response
double n_sp_taper = 6.0;            // dimensionless; use an even value such as 4, 6, or 8

inline double guiding_radius_from_L(double Lz)
{
    return Lz / v_c;
}


// Energies (inertial and rotating-frame effective)
// Not used during ramping version
double total_potential_energy(const Vec2 &cart_pos, const double &sim_time)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);

    return potential_energy_unperturbed(cart_pos) + spiral_potential(cyl_pos, sim_time, 1.0);
}

double potential_energy_unperturbed(const Vec2 &pos)
{
    return 0.5 * v_c * v_c * log(pos.x * pos.x + pos.y * pos.y);
}

// Time-dependent perturbation (ramp)
double total_PE_time(const Vec2 &cart_pos, const double &sim_time)
{
    if (pert_strength == 0.0)
        return potential_energy_unperturbed(cart_pos);

    Cyl cyl_pos = cart_to_cyl(cart_pos);

    return potential_energy_unperturbed(cart_pos)
         + spiral_potential(cyl_pos, sim_time, perturbation_scaling_factor(sim_time));
}

// Axisymmetric logarithmic potential acceleration
Vec2 LogPot_acc(const Vec2 &pos)
{
    double R2 = pos.x * pos.x + pos.y * pos.y;
    double dPhi_dx = -(v_c * v_c * pos.x) / R2;
    double dPhi_dy = -(v_c * v_c * pos.y) / R2;
    return {dPhi_dx, dPhi_dy};
}

void Leapfrog_integrator_unperturbed(Vec2 &pos, Vec2 &vel, const double &dt)
{
    // drift
    Vec2 r_half{ pos.x + 0.5 * dt * vel.x, pos.y + 0.5 * dt * vel.y };

    // kick
    Vec2 a = LogPot_acc(r_half);
    vel.x += dt * a.x;
    vel.y += dt * a.y;

    // drift
    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
}

Vec2 total_acceleration(const Vec2 &cart_pos, const double &sim_time)
{
    // If perturbation is disabled, just return axisymmetric log force
    if (pert_strength == 0.0)
        return LogPot_acc(cart_pos);

    // Guard against R -> 0.
    if (std::hypot(cart_pos.x, cart_pos.y) < 1e-12)
        return {0.0, 0.0};

    // Acceleration = -gradient(Phi).
    return LogPot_acc(cart_pos)
         - spiral_gradient(cart_pos, sim_time, perturbation_scaling_factor(sim_time));
}


void Leapfrog_integrator_perturbed(Vec2 &pos, Vec2 &vel, const double &dt, const double &sim_time)
{
    // drift
    pos.x += 0.5 * dt * vel.x;
    pos.y += 0.5 * dt * vel.y;

    // kick
    Vec2 a = total_acceleration(pos, sim_time + (0.5 * dt));
    vel.x += dt * a.x;
    vel.y += dt * a.y;

    // drift
    pos.x += 0.5 * dt * vel.x;
    pos.y += 0.5 * dt * vel.y;
}

double angular_momentum(const Vec2 &pos, const Vec2 &vel)
{
    return pos.x * vel.y - pos.y * vel.x;
}

double kinetic_energy(const Vec2 &vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y);
}

double jacobi_integral(const Vec2 &cart_pos, const Vec2 &cart_vel, const double &sim_time, const double &total_time)
{
    double KE = kinetic_energy(cart_vel);
    double PE = total_PE_time(cart_pos, sim_time);
    double Lz = angular_momentum(cart_pos, cart_vel);

    // E_J = E - Ω_p Lz  (rotating-frame constant when Φ is steady in that frame)
    return KE + PE - (omega_p * Lz);
}

void simulate_trajectory(Vec2 &orbit_pos, Vec2 &orbit_vel, const std::string &basename, double &total_time, double &dt, bool use_pert)
{
    std::ofstream trajectory_file (basename + "_position.dat");
    std::ofstream parameters_file (basename + "_parameters.dat");
    std::ofstream accel_file      (basename + "_accelerations.dat");
    std::ofstream rot_file        (basename + "_position_rot.dat");

    if (!trajectory_file || !parameters_file || !accel_file || !rot_file)
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }

    trajectory_file << "# t\tx\ty\tR\tphi\n";
    parameters_file << "# t\tL_z\tKE\tPE\tE\tE_J\tR_L\n";
    accel_file      << "# t\tax_log\tay_log\ttax_pert\tay_pert\n";
    rot_file        << "# t\tx_rot\ty_rot\n";

    int steps = static_cast<int>(total_time / dt);
    int print_interval = std::max(1, steps / 10);

    for (int i = 0; i <= steps; i++)
    {
        double sim_time = static_cast<double>(i) * dt;
        double scale = perturbation_scaling_factor(sim_time);
        bool perturb_active = use_pert && pert_strength > 0.0 && scale > 0.0;

        // Accelerations
        Vec2 log_acc = LogPot_acc(orbit_pos);
        Vec2 tot_acc = perturb_active ? total_acceleration(orbit_pos, sim_time) : log_acc;

        // Both integrators evolve inertial velocities
        double Lz = angular_momentum(orbit_pos, orbit_vel);
        double KE = kinetic_energy(orbit_vel);
        double PE = perturb_active ? total_PE_time(orbit_pos, sim_time) : potential_energy_unperturbed(orbit_pos);

        double E  = KE + PE;
        // Use the same rotating-frame diagnostic before, during, and after the ramp.
        double EJ = E - omega_p * Lz;
        double R_L = Lz / v_c;

        // Rotating frame position
        Cyl pos_cyl = cart_to_cyl(orbit_pos);
        Cyl pos_rot = pos_cyl;
        pos_rot.phi -= omega_p * sim_time;
        if (pos_rot.phi < 0)            pos_rot.phi += 2.0 * M_PI;
        if (pos_rot.phi >= 2.0 * M_PI)  pos_rot.phi -= 2.0 * M_PI;
        Vec2 pos_rot_cart = cyl_to_cart(pos_rot);

        // Output
        trajectory_file << sim_time << "\t"<< orbit_pos.x << "\t" << orbit_pos.y << "\t" << pos_cyl.R << "\t" << pos_cyl.phi << "\n";

        parameters_file << sim_time << "\t" << Lz << "\t" << KE << "\t" << PE << "\t" << E << "\t" << EJ << "\t" << R_L << "\n";

        accel_file << sim_time << "\t" << log_acc.x << "\t" << log_acc.y << "\t" << tot_acc.x << "\t" << tot_acc.y << "\n";

        rot_file << sim_time << "\t" << pos_rot_cart.x << "\t" << pos_rot_cart.y << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "[Orbit: " << basename << "] t = " << sim_time << " | L_z = " << Lz << " | E = " << E << " | E_J = " << EJ << "\n";
        }

        // The final row describes the returned state; do not advance past it.
        if (i == steps)
            break;

        if (perturb_active)
            Leapfrog_integrator_perturbed(orbit_pos, orbit_vel, dt, sim_time);
        else
            Leapfrog_integrator_unperturbed(orbit_pos, orbit_vel, dt);
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
    rot_file.close();
}

// Compare force due to logarithmic potential and spiral perturbation
void map_force_ratios(double xmin, double xmax, double ymin, double ymax, double step, double sim_time, const double &total_time, const std::string& filename)
{
    std::ofstream fout(filename);
    fout << "# x\ty\t|F_perturb|/|F_axisym|\n";

    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos{ x, y };

            Vec2 F_axisym = LogPot_acc(pos);
            Vec2 F_total  = total_acceleration(pos, sim_time);

            Vec2 F_perturb{ F_total.x - F_axisym.x,  F_total.y - F_axisym.y };

            double mag_axisym = std::hypot(F_axisym.x, F_axisym.y);
            double mag_pert   = std::hypot(F_perturb.x, F_perturb.y);

            double ratio = (mag_axisym > 1e-12) ? (mag_pert / mag_axisym) : 0.0;

            fout << x << "\t" << y << "\t" << ratio << "\n";
        }
        fout << "\n";
    }

    fout.close();
    std::cout << "Wrote force ratio map to: " << filename << "\n";
}

// Map spiral potential to total potential
void map_potential_ratios(double xmin, double xmax, double ymin, double ymax, double step, double sim_time, const double& total_time, const std::string& filename)
{
    std::ofstream fout(filename);
    fout << "# x\ty\t|Phi_pert|/|Phi_tot|\n";

    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos{ x, y };
            Cyl cyl_pos = cart_to_cyl(pos);

            double phi_log = 0.5 * v_c * v_c * std::log(pos.x * pos.x + pos.y * pos.y);

            double phi_pert = spiral_potential(cyl_pos, sim_time, perturbation_scaling_factor(sim_time));

            double phi_tot = phi_log + phi_pert;
            double ratio = (std::abs(phi_tot) > 1e-12) ? std::abs(phi_pert) / std::abs(phi_tot) : 0.0;

            fout << x << '\t' << y << '\t' << ratio << '\n';
        }
        fout << '\n';
    }

    fout.close();
    std::cout << "Wrote potential ratio map to: " << filename << "\n";
}

// JR and random energy
double getJRunperturbed(Vec2 pos0, Vec2 vel, double dt)
{
    double JRback = 0.0;

    // Initial cylindrical coordinates
    double Rold = std::hypot(pos0.x, pos0.y);
    double vold = (pos0.x * vel.x + pos0.y * vel.y) / Rold; // v_R

    // First leapfrog step
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    double Rmid = std::hypot(pos0.x, pos0.y);
    double vmid = (pos0.x * vel.x + pos0.y * vel.y) / Rmid;

    // Second leapfrog step
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    double Rnew = std::hypot(pos0.x, pos0.y);
    double vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    // Find first extremum (sign change of discrete second derivative)
    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;
    }

    // move away from extremum
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    Rold = Rmid; vold = vmid;
    Rmid = Rnew; vmid = vnew;
    Rnew = std::hypot(pos0.x, pos0.y);
    vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    // next extremum
    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

        JRback += (Rnew - Rmid) * 0.5 * (vnew + vmid);
    }

    // Continue to next extremum
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    Rold = Rmid; vold = vmid;
    Rmid = Rnew; vmid = vnew;
    Rnew = std::hypot(pos0.x, pos0.y);
    vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    JRback += (Rnew - Rmid) * 0.5 * (vnew + vmid);

    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

        JRback += (Rnew - Rmid) * 0.5 * (vnew + vmid);
    }

    return JRback / (2.0 * M_PI); // J_R
}

double random_energy(const Vec2 &pos, const Vec2 &vel)
{
    double Lz = angular_momentum(pos, vel);
    double Rg = Lz / v_c;

    double R = std::hypot(pos.x, pos.y);
    double v_R = (pos.x * vel.x + pos.y * vel.y) / R;

    double potential_R  = 0.5 * v_c * v_c * log(R  * R);
    double potential_Rg = 0.5 * v_c * v_c * log(Rg * Rg);

    double eff_pot_R  = potential_R  + (Lz * Lz) / (2.0 * R  * R);
    double eff_pot_Rg = potential_Rg + (Lz * Lz) / (2.0 * Rg * Rg);

    double eR = 0.5 * v_R * v_R + eff_pot_R - eff_pot_Rg;
    return eR;
}

void integrate_no_io(Vec2 pos, Vec2 vel, double total_time, double dt, bool use_pert, double& Lz0, double& Lz_end,
                            double& JR0, double& JR_end, double& eR0, double& eR_end, Vec2* pos_end, Vec2* vel_end)
{
    Lz0 = angular_momentum(pos, vel);
    JR0 = getJRunperturbed(pos, vel, dt);
    eR0 = random_energy(pos, vel);

    int steps = static_cast<int>(total_time / dt);
    for (int i = 0; i < steps; i++)
    {
        const double t = i * dt;
        const bool perturb_active = use_pert && pert_strength > 0.0 && (perturbation_scaling_factor(t) > 0.0);

        if (perturb_active)
            Leapfrog_integrator_perturbed(pos, vel, dt, t);
        else
            Leapfrog_integrator_unperturbed(pos, vel, dt);
    }

    Lz_end = angular_momentum(pos, vel);
    JR_end = getJRunperturbed(pos, vel, dt);
    eR_end = random_energy(pos, vel);

    if (pos_end) *pos_end = pos;
    if (vel_end) *vel_end = vel;
}

/*
// Finite-difference acceleration from total_PE_time
Vec2 total_acceleration(const Vec2 &cart_pos, const double &sim_time, const double &total_time)
{
    if (pert_strength == 0.0)
    {
        return LogPot_acc(cart_pos);
    }
    double delta = 1e-6;

    Vec2 pos_dx_p{ cart_pos.x + delta, cart_pos.y };
    Vec2 pos_dx_m{ cart_pos.x - delta, cart_pos.y };
    Vec2 pos_dy_p{ cart_pos.x, cart_pos.y + delta };
    Vec2 pos_dy_m{ cart_pos.x, cart_pos.y - delta };

    double dVdx = (total_PE_time(pos_dx_m, sim_time) - total_PE_time(pos_dx_p, sim_time)) * 0.5 / delta;
    double dVdy = (total_PE_time(pos_dy_m, sim_time) - total_PE_time(pos_dy_p, sim_time)) * 0.5 / delta;

    return { dVdx, dVdy };
}
*/


  
