#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <functional>
#include "VecUtils.h"
#include "LogPot.h"

using std::string;
using std::vector;
using std::ofstream;
using GradCandidate = std::tuple<double, double, double>; // [x, y, gradient]
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


double spiral_radial_envelope(double R)
{
    if (R_sp_taper <= 0.0)
        return 1.0;

    if (R <= 0.0)
        return 0.0;

    const double q  = R / R_sp_taper;
    const double qn = std::pow(q, n_sp_taper);
    return qn / (1.0 + qn);
}

double d_spiral_radial_envelope_dR(double R)
{
    if (R_sp_taper <= 0.0 || R <= 1e-12)
        return 0.0;

    const double W = spiral_radial_envelope(R);
    return (n_sp_taper / R) * W * (1.0 - W);
}

// Effective potential 
double effective_potential(const Vec2 &cart_pos, const bool &use_pert)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);

    // Axisymmetric log potential
    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    // Optional spiral perturbation
    double pert_energy = 0.0;
    if (use_pert)
    {
        double kappa = alpha / cyl_pos.R;
        double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
        double W = spiral_radial_envelope(cyl_pos.R);

        // constant Spiral-phase
        pert_energy = (Phi_S * W * pert_strength) * cos((alpha * log(cyl_pos.R / R_CR)) - m * cyl_pos.phi);

        // For a bar:
        // pert_energy = (Phi_S * pert_strength) * cos(m * cyl_pos.phi);
    }

    // Phi_eff = Phi - 0.5 * omega_p^2 * R^2 (Binney & Tremaine, eq. 3.114).
    double cent_term = 0.5 * omega_p * omega_p * (cyl_pos.R * cyl_pos.R);

    return log_pot_energy + pert_energy - cent_term;
}

Vec2 effective_potential_gradient(const Vec2 &cart_pos, const bool &use_pert)
{
    double delta = 1e-8;

    Vec2 dx_pos = {cart_pos.x + delta, cart_pos.y};
    Vec2 dx_neg = {cart_pos.x - delta, cart_pos.y};
    Vec2 dy_pos = {cart_pos.x, cart_pos.y + delta};
    Vec2 dy_neg = {cart_pos.x, cart_pos.y - delta};

    Vec2 grad;
    grad.x = (effective_potential(dx_pos, use_pert) - effective_potential(dx_neg, use_pert)) / (2.0 * delta);
    grad.y = (effective_potential(dy_pos, use_pert) - effective_potential(dy_neg, use_pert)) / (2.0 * delta);

    return grad;
}

// Lagrange point identification
bool is_near_cluster(const vector<GradCandidate>& cluster, const GradCandidate& point, double epsilon)
{
    for (const auto& c : cluster)
    {
        double dx = std::get<0>(c) - std::get<0>(point);
        double dy = std::get<1>(c) - std::get<1>(point);
        if (std::sqrt(dx * dx + dy * dy) < epsilon) return true;
    }
    return false;
}

void scan_Lagrange_points(double xmin, double xmax, double ymin, double ymax, double step, const string& filename, bool use_pert)
{
    std::ofstream fout(filename);
    std::ofstream fout_lagrange("lagrange_points.dat");

    fout << "# x\ty\t|∇Phi_eff|\n";
    fout_lagrange << "# Lagrange points\n";

    vector<GradCandidate> low_grad_points;
    double threshold = 5e-5;
    double epsilon   = 0.5; // clustering radius in kpc

    // Find all gradients below the threshold
    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos = {x, y};
            Vec2 grad = effective_potential_gradient(pos, use_pert);
            double grad_mag = std::sqrt(grad.x * grad.x + grad.y * grad.y);

            fout << x << " " << y << " " << grad_mag << "\n";

            if (grad_mag < threshold)
                low_grad_points.emplace_back(x, y, grad_mag);
        }
        fout << "\n";
    }

    // Sort by gradient magnitude
    std::sort(low_grad_points.begin(), low_grad_points.end(), [](const GradCandidate& a, const GradCandidate& b)
              { return std::get<2>(a) < std::get<2>(b); });

    // Cluster
    vector<vector<GradCandidate>> clusters;
    for (const auto& pt : low_grad_points)
    {
        bool added = false;
        for (auto& cluster : clusters)
        {
            if (is_near_cluster(cluster, pt, epsilon))
            {
                cluster.push_back(pt);
                added = true;
                break;
            }
        }
        if (!added)
            clusters.push_back({pt});
    }

    // Output one representative (min grad) per cluster
    for (const auto& cluster : clusters)
    {
        auto min_point = *std::min_element(cluster.begin(), cluster.end(), [](const GradCandidate& a, const GradCandidate& b)
            { return std::get<2>(a) < std::get<2>(b); });

        fout_lagrange << std::get<0>(min_point) << " " << std::get<1>(min_point) << " " << std::get<2>(min_point) << "\n";
    }

    fout.close();
    fout_lagrange.close();
}

double perturbation_scaling_factor(double sim_time)
{
    //Variables for perturbation profiless
    double t_start  = 2000.0;
    double dur_up   = 500.0;
    double dur_hold = 0.0;
    double dur_down = 500.0;
    double t_up_end   = t_start + dur_up;
    double t_hold_end = t_up_end + dur_hold;
    double t_down_end = t_hold_end + dur_down;

    if (sim_time < t_start)
        return 0.0;

    // Ramp up
    if (dur_up > 0.0 && sim_time < t_up_end)
        return (sim_time - t_start) / dur_up;

    // Immediate drop if dur_up == 0
    if (dur_up == 0.0 && sim_time < t_hold_end)
        return 1.0;

    // Hold
    if (sim_time < t_hold_end)
        return 1.0;

    // Ramp down
    if (dur_down > 0.0 && sim_time < t_down_end)
        return (t_down_end - sim_time) / dur_down;
    
    return 0.0;
}

// Energies (inertial and rotating-frame effective)
// Not used during ramping version
double total_potential_energy(const Vec2 &cart_pos, const double &sim_time)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);

    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
    double W = spiral_radial_envelope(cyl_pos.R);

    double pert_energy = (Phi_S * W * pert_strength) * cos(alpha * log(cyl_pos.R / R_CR) + m * omega_p * sim_time - m * cyl_pos.phi);

    // For a bar:
    // double pert_energy = (Phi_S * pert_strength) * cos(m*omega_p*sim_time + m * cyl_pos.phi);

    return log_pot_energy + pert_energy;
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

    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
    double W = spiral_radial_envelope(cyl_pos.R);

    double scale = perturbation_scaling_factor(sim_time);

    double pert_energy = (Phi_S * W * pert_strength * scale) * cos(alpha * log(cyl_pos.R / R_CR) + m * omega_p * sim_time - m * cyl_pos.phi);

    // For a bar:
    // double pert_energy = (Phi_S * pert_strength * scale) * cos(m * omega_p * sim_time - m * cyl_pos.phi);

    return log_pot_energy + pert_energy;
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

    double x = cart_pos.x;
    double y = cart_pos.y;

    double R_squared = x*x + y*y;
    double R  = std::sqrt(R_squared);
    double phi = std::atan2(y, x);

    //Guard against R -> 0 
    if (R < 1e-12)
        return {0.0, 0.0};

    double dPhi0_dR = (v_c * v_c) / R;

    double Sigma = Sigma_o * std::exp(-R / R_o);
    double kappa = alpha / R;
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;

    //Find the amplitude(A) of the perturbation
    double scale = perturbation_scaling_factor(sim_time);

    double W    = spiral_radial_envelope(R);
    double dW_dR = d_spiral_radial_envelope_dR(R);

    double pert_Amplitude = Phi_S * W * pert_strength * scale;

    //ψ = alpha ln(R/R_CR) + m ω_p t - mφ
    double psi = alpha * std::log(R / R_CR) + m * omega_p * sim_time - m * phi;

    double cpsi = cos(psi);
    double spsi = sin(psi);

    // d/dR [Phi_S(R) * W(R) * pert_strength * scale]
    // Phi_S ∝ R exp(-R/R_o), so dPhi_S/dR = Phi_S * (1/R - 1/R_o).
    double dPhiS_dR = Phi_S * (1.0/R - 1.0/R_o);
    double dA_dR = pert_strength * scale * (dPhiS_dR * W + Phi_S * dW_dR);

    double dpsi_dR   = alpha / R;
    double dpsi_dphi = -m;

    double dPhi_dR   = dPhi0_dR + dA_dR * cpsi - pert_Amplitude * spsi * dpsi_dR;
    double dPhi_dphi = m * pert_Amplitude * spsi; 

    double dPhi_dx = dPhi_dR * (x / R) + dPhi_dphi * (-y / R_squared);
    double dPhi_dy = dPhi_dR * (y / R) + dPhi_dphi * ( x / R_squared);

    //Acceleration=-∇Φ
    return { -dPhi_dx, -dPhi_dy };
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
    if (pert_strength == 0.0)
    {
        double KE = kinetic_energy(cart_vel);
        double PE = potential_energy_unperturbed(cart_pos);
        return KE + PE;
    }

    double KE = kinetic_energy(cart_vel);
    double PE = total_PE_time(cart_pos, sim_time);
    double Lz = angular_momentum(cart_pos, cart_vel);

    // E_J = E - Ω_p Lz  (rotating-frame constant when Φ is steady in that frame)
    return KE + PE - (omega_p * Lz);
}

double effective_potential_time(const Vec2 &cart_pos, const bool &use_pert, double sim_time)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);

    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double pert_energy = 0.0;
    if (use_pert && pert_strength != 0.0)
    {
        double kappa = alpha / cyl_pos.R;
        double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
        double W = spiral_radial_envelope(cyl_pos.R);
        double scale = perturbation_scaling_factor(sim_time);

        pert_energy = (Phi_S * W * pert_strength * scale) * cos((alpha * log(cyl_pos.R / R_CR)) + m * omega_p * sim_time - m * cyl_pos.phi);
        // For a bar:
        // pert_energy = (Phi_S * pert_strength * scale) * cos(m * omega_p * sim_time - m * cyl_pos.phi);
    }

    double cent_term = 0.5 * omega_p * omega_p * (cyl_pos.R * cyl_pos.R);

    return log_pot_energy + pert_energy - cent_term;
}

Vec2 effective_potential_gradient_time(const Vec2 &cart_pos, const bool &use_pert, double sim_time, double total_time)
{
    double delta = 1e-8;

    Vec2 dx_pos = {cart_pos.x + delta, cart_pos.y};
    Vec2 dx_neg = {cart_pos.x - delta, cart_pos.y};
    Vec2 dy_pos = {cart_pos.x, cart_pos.y + delta};
    Vec2 dy_neg = {cart_pos.x, cart_pos.y - delta};

    Vec2 grad;
    grad.x = (effective_potential_time(dx_pos, use_pert, sim_time) - effective_potential_time(dx_neg, use_pert, sim_time)) / (2.0 * delta);
    grad.y = (effective_potential_time(dy_pos, use_pert, sim_time) -  effective_potential_time(dy_neg, use_pert, sim_time)) / (2.0 * delta);
              
    return grad;
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
    int print_interval = steps / 10;

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
        double EJ = perturb_active ? jacobi_integral(orbit_pos, orbit_vel, sim_time, total_time) : E;
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

            double kappa = alpha / cyl_pos.R;
            double Sigma = Sigma_o * std::exp(-cyl_pos.R / R_o);
            double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
            double W = spiral_radial_envelope(cyl_pos.R);

            double scale = perturbation_scaling_factor(sim_time);

            double phi_pert = (Phi_S * W * pert_strength * scale) * std::cos(alpha * std::log(cyl_pos.R / R_CR) + m * omega_p * sim_time - m * cyl_pos.phi);

            // For a bar:
            // double phi_pert = (Phi_S * pert_strength * scale) * std::cos(m * cyl_pos.phi);

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

static void integrate_no_io(Vec2 pos, Vec2 vel, double total_time, double dt, bool use_pert, double& Lz0, double& Lz_end, 
                            double& JR0, double& JR_end, double& eR0, double& eR_end, Vec2* pos_end = nullptr, Vec2* vel_end = nullptr)
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

int main(int argc, char** argv)
{
    bool use_pert = false;
    if (argc > 1)
    {
        std::string arg = argv[1];
        if      (arg == "-pert")     use_pert = true;
        else if (arg == "-no_pert")  use_pert = false;
        else
        {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        }
    }

    // Simulation durations
    double total_time = 5000.0;
    double dt = 0.01;

    // --- Time series diagnostics at R=10 kpc, phi=0 ---
    std::ofstream fpert("phi_perturbation_time_series.dat");
    std::ofstream frat ("phi_ratio_time_series.dat");
    std::ofstream fstr ("pert_strength_time_series.dat");

    if (!fpert || !frat || !fstr)
    {
        std::cerr << "Error: could not open one of the time-series output files.\n";
        return 1;
    }

    fpert << "# time [Myr]\tPhi_perturb\n";
    frat  << "# time [Myr]\t|Phi_pert| / |Phi_total|\n";
    fstr  << "# time [Myr]\ts(t)=ramp(0..1)\tA(t)=pert_strength*s(t)\n";

    Vec2 pos_probe = cyl_to_cart({10.0, 0.0});
    for (double t = 0.0; t <= total_time; t += dt)
    {
        Cyl cyl = cart_to_cyl(pos_probe);

        // Guard (shouldn't happen here, but keeps formulas safe)
        if (cyl.R < 1e-12)
        {
            fpert << t << "\t0\n";
            frat  << t << "\t0\n";
            fstr  << t << "\t0\t0\n";
            continue;
        }

        double kappa = alpha / cyl.R;
        double Sigma = Sigma_o * std::exp(-cyl.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
        double W = spiral_radial_envelope(cyl.R);

        double scale = perturbation_scaling_factor(t);

        double phi_pert  = (Phi_S * W * pert_strength * scale) * cos(alpha * log(cyl.R / R_CR) + m * omega_p * t - m * cyl.phi);

        double phi_total = 0.5 * v_c * v_c * std::log(cyl.R * cyl.R) + phi_pert;

        fpert << t << "\t" << phi_pert << "\n";

        double ratio = (std::fabs(phi_total) > 1e-12) ? (std::fabs(phi_pert) / std::fabs(phi_total)) : 0.0;
        frat  << t << "\t" << ratio << "\n";

        fstr  << t << "\t" << scale << "\t" << (pert_strength * scale) << "\n";
    }

    fpert.close();
    frat.close();
    fstr.close();

    std::cout << "Wrote phi_perturbation_time_series.dat, phi_ratio_time_series.dat, and pert_strength_time_series.dat\n";

    double x_min = -20.0, x_max = 20.0;
    double y_min = -20.0, y_max = 20.0;
    double step  = 0.1;

    {
        std::ofstream fout("phi_eff.dat");
        if (!fout)
        {
            std::cerr << "Error: could not open phi_eff.dat\n";
            return 1;
        }

        fout << "# x [kpc]\ty [kpc]\tPhi_eff(x,y)\n";

        const double snapshot_time = 2500.0; // peak of ramp
        for (double y = y_min; y <= y_max; y += step)
        {
            for (double x = x_min; x <= x_max; x += step)
            {
                Vec2 pos_xy{ x, y };
                double phi_eff = effective_potential_time(pos_xy, use_pert, snapshot_time);
                fout << x << "\t" << y << "\t" << phi_eff << "\n";
            }
            fout << "\n";
        }
    }

    std::cout << "Wrote phi_eff.dat (perturbation " << (use_pert ? "on" : "off") << ")\n";

    scan_Lagrange_points(x_min, x_max, y_min, y_max, step, "grad_phi_eff.dat", use_pert);
    std::cout << "Wrote grad_phi_eff.dat and lagrange_points.dat\n";

    {
        std::ifstream fin("DF_initial_conditions.dat");
        if (!fin)
        {
            std::cerr << "Could not open DF_initial_conditions.dat\n";
        }
        else
        {
            std::ofstream out("deltaLz_from_DF.dat");
            if (!out)
            {
                std::cerr << "Could not open deltaLz_from_DF.dat for writing\n";
                return 1;
            }

            // DF file format (6 columns):
            // R  L  v_R  v_phi  x  y
            out << "# id  R  x  y  vR  vphi  Lz0  Lz_end  dLz  JR0  JR_end  dJR  eR0  eR_end  deR  vR_end\n";

            std::string line;
            int id = 0;
            std::size_t N_rows = 0;
            double sum_dLz = 0.0, sum_dJR = 0.0;

            while (std::getline(fin, line))
            {
                if (line.empty() || line[0] == '#') continue;

                double R, L, vR, vphi, x, y;
                std::istringstream iss(line);
                if (!(iss >> R >> L >> vR >> vphi >> x >> y))
                    continue;

                Vec2 pos0{ x, y };
                Vec2 vel0 = cyl_to_cart_vel(vR, vphi, pos0);

                Vec2 pos_end, vel_end;
                double Lz0 = 0.0, Lz_end = 0.0;
                double JR0 = 0.0, JR_end = 0.0;
                double eR0 = 0.0, eR_end = 0.0;

                integrate_no_io(pos0, vel0, total_time, dt, use_pert, Lz0, Lz_end, JR0, JR_end, eR0, eR_end, &pos_end, &vel_end);

                const double dLz = (Lz_end - Lz0);
                const double dJR = (JR_end - JR0);
                const double deR = (eR_end - eR0);

                const double vR_end = cart_to_cyl_inertial(pos_end, vel_end).x;

                out << id << ' ' << R << ' ' << x << ' ' << y << ' ' << vR << ' ' << vphi << ' ' << Lz0 << ' ' << Lz_end << ' ' << dLz << ' '
                    << JR0 << ' ' << JR_end << ' ' << dJR << ' ' << eR0 << ' ' << eR_end << ' ' << deR << ' ' << vR_end << '\n';

                sum_dLz += dLz;
                sum_dJR += dJR;
                id++;
                N_rows++;
            }

            out.close();
            fin.close();

            double mean_dLz = (N_rows > 0) ? (sum_dLz / N_rows) : 0.0;
            std::cout << "Average ΔLz = " << mean_dLz << " (from " << N_rows << " stars)\n";
            std::cout << "Wrote deltaLz_from_DF.dat (6-col DF input)\n";
        }
    }

    {
        std::ifstream fin2("deltaLz_from_DF.dat");
        if (!fin2)
        {
            std::cerr << "Could not open deltaLz_from_DF.dat for RMS calc\n";
        }
        else
        {
            double sum2 = 0.0;
            int N = 0;

            std::string line;
            while (std::getline(fin2, line))
            {
                if (line.empty() || line[0] == '#') continue;

                int id;
                double R, x, y, vR, vphi;
                double Lz0, Lz_end, dLz;
                double JR0, JR_end, dJR;
                double eR0, eR_end, deR;
                double vR_end;

                std::istringstream iss(line);
                if (!(iss >> id >> R >> x >> y >> vR >> vphi
                          >> Lz0 >> Lz_end >> dLz
                          >> JR0 >> JR_end >> dJR
                          >> eR0 >> eR_end >> deR
                          >> vR_end))
                    continue;

                if (Lz0 >= 2.0 && Lz0 <= 16.0) 
                {
                    sum2 += dLz * dLz;
                    ++N;
                }
            }

            if (N > 0)
            {
                double rms = std::sqrt(sum2 / double(N));
                std::cout << "RMS(ΔLz) in CR window [2.0, 16.0] = " << rms << " (N=" << N << ")\n";
            }
            else
            {
                std::cout << "No stars in CR window [2.0, 16.0]\n";
            }
        }
    }

    //Single orbits for parameter verification
    std::vector<Vec2> initial_pos =
    {
        { 10.2,  0.0 },
        { -8.8, -4.5 },
        { -4.6,  8.8 },
        {  8.6,  4.4 }
    };

    for (size_t i = 0; i < initial_pos.size(); i++)
    {
        Vec2 p = initial_pos[i];
        double phi = std::atan2(p.y, p.x);

        double v_phi = 0.98 * v_c;
        double v_R   = 0.02 * v_c;

        Vec2 v;
        v.x = v_R * std::cos(phi) - v_phi * std::sin(phi);
        v.y = v_R * std::sin(phi) + v_phi * std::cos(phi);

        std::ostringstream basename;
        basename << "orbit_L" << i;
        simulate_trajectory(p, v, basename.str(), total_time, dt, use_pert);
    }
    return 0;
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


  
