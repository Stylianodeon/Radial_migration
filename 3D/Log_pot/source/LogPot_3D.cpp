#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "LogPot_3D.h"

using std::string;
using std::vector;

const double PI = std::acos(-1.0);

// Global axisymmetric potential parameters
double v_c = 0.22;     // kpc/Myr = 220 km/s
double q   = 0.9;      // flattening parameter

// Global spiral perturbation parameters
double m           = 4.0;           // spiral mode
double G           = 4.5e-12;       // kpc^3 / (Msun Myr^2)
double R_o         = 2.5;           // kpc
double R_CR        = 10.0;          // kpc
double Sigma_o     = 1.226e+9;      // Msun / kpc^2
double e_s         = 0.1;           // fractional surface-density spiral amplitude
double theta_deg   = 15.0;          // pitch angle in degrees
double theta       = theta_deg * PI / 180.0;
double alpha       = m / std::tan(theta);
double omega_p     = 0.022;         // pattern speed in rad/Myr
double pert_strength = 0.1;         // additional numerical scaling factor

// Vertical scale height. Large values make the perturbation almost independent of z.
double spiral_z_scale = 0.3;        // kpc

// Time envelope for the perturbation.
double t_spiral_start     = 2000.0; // Myr
double t_spiral_ramp_up   = 500.0;  // Myr
double t_spiral_hold      = 0.0;    // Myr
double t_spiral_ramp_down = 500.0;  // Myr


// Phi_0 = 0.5 * v_c^2 * ln(x^2 + y^2 + z^2/q^2)
double potential_energy(const Vec3& pos)
{
    const double q2 = q * q;
    const double S = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (S <= 0.0)
    {
        std::cerr << "Error: logarithmic potential is singular at the origin.\n";
        return 0.0;
    }

    return 0.5 * v_c * v_c * std::log(S);
}

Vec3 LogPot_acc(const Vec3& pos)
{
    const double q2 = q * q;
    const double S = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (S <= 0.0)
    {
        std::cerr << "Error: acceleration is singular at the origin.\n";
        return {0.0, 0.0, 0.0};
    }

    const double ax = -v_c * v_c * pos.x / S;
    const double ay = -v_c * v_c * pos.y / S;
    const double az = -v_c * v_c * pos.z / (q2 * S);

    return {ax, ay, az};
}

// ======================================================
// Spiral perturbation
// Phi_sp(R,phi,z,t) = A(R,t) Z(z) cos(psi)
//
// A(R,t)  = pert_strength * s(t) * Phi_S(R)
// Phi_S   = 2 pi G e_s Sigma(R) / |k(R)|
// Sigma   = Sigma_o exp(-R/R_o)
// k(R)    = alpha / R
// psi     = alpha ln(R/R_CR) + m omega_p t - m phi
// Z(z)    = sech^2(z / spiral_z_scale)
// s(t)    = perturbation scaling factor / time envelope
// ======================================================

double perturbation_scaling_factor(double sim_time)
{
    const double t_up_end   = t_spiral_start + t_spiral_ramp_up;
    const double t_hold_end = t_up_end + t_spiral_hold;
    const double t_down_end = t_hold_end + t_spiral_ramp_down;

    if (sim_time < t_spiral_start)
        return 0.0;

    if (t_spiral_ramp_up > 0.0 && sim_time < t_up_end)
        return (sim_time - t_spiral_start) / t_spiral_ramp_up;

    if (t_spiral_ramp_up == 0.0 && sim_time < t_hold_end)
        return 1.0;

    if (sim_time < t_hold_end)
        return 1.0;

    if (t_spiral_ramp_down > 0.0 && sim_time < t_down_end)
        return (t_down_end - sim_time) / t_spiral_ramp_down;

    return 0.0;
}

// Smooth vertical attenuation of the spiral perturbation.
// The spiral is strongest in the disc midplane, Z(0)=1, and weakens
// away from the plane. Setting spiral_z_scale very large gives Z~1
// and recovers an effectively 2D perturbation.
double spiral_vertical_envelope(double z)
{
    if (spiral_z_scale <= 0.0)
        return 1.0;

    const double u = z / spiral_z_scale;

    // Avoid overflow in cosh for very large |z|.
    if (std::abs(u) > 20.0)
        return 0.0;

    // sech^2(u) = 1 / cosh^2(u).
    const double c = std::cosh(u);
    return 1.0 / (c * c);
}

double d_spiral_vertical_envelope_dz(double z)
{
    if (spiral_z_scale <= 0.0)
        return 0.0;

    const double Z = spiral_vertical_envelope(z);
    const double u = z / spiral_z_scale;

    // d/dz sech^2(z/h_z) = -(2/h_z) sech^2(z/h_z) tanh(z/h_z).
    return -(2.0 / spiral_z_scale) * Z * std::tanh(u);
}

double spiral_potential_3D(const Vec3& pos, double sim_time)
{
    const double x = pos.x, y = pos.y, z = pos.z;

    const double R2 = x * x + y * y;
    const double R  = std::sqrt(R2);

    if (R < 1e-12)
        return 0.0;

    const double phi = std::atan2(y, x);

    const double Sigma = Sigma_o * std::exp(-R / R_o);
    const double k_abs = std::abs(alpha) / R;
    const double Phi_S = (2.0 * PI * G * e_s * Sigma) / k_abs;

    const double Z     = spiral_vertical_envelope(z);
    const double scale = perturbation_scaling_factor(sim_time);

    const double A = Phi_S * pert_strength * scale;

    const double psi = alpha * std::log(R / R_CR) + m * omega_p * sim_time - m * phi;

    return A * Z * std::cos(psi);
}

Vec3 spiral_acceleration_3D(const Vec3& pos, double sim_time)
{
    const double x = pos.x, y = pos.y, z = pos.z;

    const double R2 = x * x + y * y;
    const double R  = std::sqrt(R2);

    if (R < 1e-12 || pert_strength == 0.0)
        return {0.0, 0.0, 0.0};

    const double phi = std::atan2(y, x);

    const double Sigma = Sigma_o * std::exp(-R / R_o);
    const double k_abs = std::abs(alpha) / R;
    const double Phi_S = (2.0 * PI * G * e_s * Sigma) / k_abs;

    const double scale = perturbation_scaling_factor(sim_time);
    if (scale == 0.0)
        return {0.0, 0.0, 0.0};

    const double Z     = spiral_vertical_envelope(z);
    const double dZ_dz = d_spiral_vertical_envelope_dz(z);

    // A(R,t) excludes the vertical envelope.
    const double A = Phi_S * pert_strength * scale;

    // Phi_S proportional to R exp(-R/R_o), so
    // dPhi_S/dR = Phi_S * (1/R - 1/R_o).
    const double dPhiS_dR = Phi_S * (1.0 / R - 1.0 / R_o);
    const double dA_dR = pert_strength * scale * dPhiS_dR;

    const double psi = alpha * std::log(R / R_CR) + m * omega_p * sim_time - m * phi;

    const double cpsi = std::cos(psi);
    const double spsi = std::sin(psi);

    const double dpsi_dR = alpha / R;

    // Cylindrical derivatives of Phi_sp.
    const double dPhi_dR   = Z * (dA_dR * cpsi - A * spsi * dpsi_dR);
    const double dPhi_dphi = m * A * Z * spsi;
    const double dPhi_dz   = A * dZ_dz * cpsi;

    // Convert gradient from cylindrical to Cartesian.
    const double dPhi_dx = dPhi_dR * (x / R) + dPhi_dphi * (-y / R2);
    const double dPhi_dy = dPhi_dR * (y / R) + dPhi_dphi * ( x / R2);

    return {-dPhi_dx, -dPhi_dy, -dPhi_dz};
}

double total_potential_energy_3D(const Vec3& pos, double sim_time, bool use_pert)
{
    double Phi = potential_energy(pos);

    if (use_pert)
        Phi += spiral_potential_3D(pos, sim_time);

    return Phi;
}

Vec3 total_acceleration_3D(const Vec3& pos, double sim_time, bool use_pert)
{
    Vec3 a = LogPot_acc(pos);

    if (!use_pert)
        return a;

    Vec3 asp = spiral_acceleration_3D(pos, sim_time);

    return {a.x + asp.x, a.y + asp.y, a.z + asp.z};
}

// ======================================================
// Initial-condition file reader
// ======================================================
struct InitialStar3D
{
    int id;

    double R;
    double phi;
    double z;

    double v_R;
    double v_phi;
    double v_z;

    double Lz;
    double Rg;

    Vec3 pos;
    Vec3 vel;
};

std::vector<InitialStar3D> read_initial_conditions_3D(const std::string& filename, int max_stars)
{
    std::ifstream in(filename);

    std::vector<InitialStar3D> stars;
    std::string line;

    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        InitialStar3D star;

        double x, y, z_cart;
        double vx, vy, vz_cart;

        iss >> star.id >> star.R >> star.phi >> star.z >> star.v_R >> star.v_phi >> star.v_z 
            >> star.Lz >> star.Rg >> x >> y >> z_cart >> vx >> vy >> vz_cart;

        star.pos = {x, y, z_cart};
        star.vel = {vx, vy, vz_cart};

        stars.push_back(star);
    }
    return stars;
}


// ======================================================
// Integrators
// ======================================================
void Leapfrog_integrator_unperturbed(Vec3& pos, Vec3& vel, const double& dt)
{
    Vec3 r_half
    {
        pos.x + 0.5 * dt * vel.x,
        pos.y + 0.5 * dt * vel.y,
        pos.z + 0.5 * dt * vel.z
    };

    Vec3 a = LogPot_acc(r_half);

    vel.x += dt * a.x;
    vel.y += dt * a.y;
    vel.z += dt * a.z;

    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
    pos.z = r_half.z + 0.5 * dt * vel.z;
}

void Leapfrog_integrator_3D(Vec3& pos, Vec3& vel, const double& dt, double sim_time, bool use_pert)
{
    Vec3 r_half
    {
        pos.x + 0.5 * dt * vel.x,
        pos.y + 0.5 * dt * vel.y,
        pos.z + 0.5 * dt * vel.z
    };

    // Drift to the midpoint, then evaluate the force at both the midpoint
    // position and midpoint time. This is important because the spiral
    // perturbation rotates explicitly with time.
    Vec3 a = total_acceleration_3D(r_half, sim_time + 0.5 * dt, use_pert);

    vel.x += dt * a.x;
    vel.y += dt * a.y;
    vel.z += dt * a.z;

    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
    pos.z = r_half.z + 0.5 * dt * vel.z;
}

void leapfrog_step(Vec3& pos, Vec3& vel, double dt)
{
    Leapfrog_integrator_unperturbed(pos, vel, dt);
}


// ======================================================
// Diagnostics
// ======================================================
double kinetic_energy(const Vec3& vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
}

double total_energy(const Vec3& pos, const Vec3& vel)
{
    // Background/unperturbed energy only.
    return kinetic_energy(vel) + potential_energy(pos);
}

double total_energy_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert)
{
    // Total energy including the spiral perturbation if use_pert = true.
    return kinetic_energy(vel) + total_potential_energy_3D(pos, sim_time, use_pert);
}

Vec3 angular_momentum(const Vec3& pos, const Vec3& vel)
{
    return cross_product(pos, vel);
}

double angular_momentum_z(const Vec3& pos, const Vec3& vel)
{
    return pos.x * vel.y - pos.y * vel.x;
}

double jacobi_integral_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert)
{
    const double E  = total_energy_3D(pos, vel, sim_time, use_pert);
    const double Lz = angular_momentum_z(pos, vel);

    return E - omega_p * Lz;
}

double random_energy(const Vec3& pos, const Vec3& vel)
{
    // E_rand = E - E_c(Lz)
    const double E  = total_energy(pos, vel);
    const double Lz = angular_momentum_z(pos, vel);

    const double Rg = std::abs(Lz) / v_c;

    if (Rg < 1e-12)
        return 0.0;

    const double Phi_Rg = 0.5 * v_c * v_c * std::log(Rg * Rg);
    const double Ec     = Phi_Rg + 0.5 * v_c * v_c;

    return E - Ec;
}

double vertical_energy(const Vec3& pos, const Vec3& vel)
{
    // Ez = 0.5 vz^2 + Phi0(R,z) - Phi0(R,0).
    const double R = std::hypot(pos.x, pos.y);

    if (R < 1e-12)
        return 0.0;

    const Vec3 midplane_pos{R, 0.0, 0.0};

    const double Phi_Rz = potential_energy(pos);
    const double Phi_R0 = potential_energy(midplane_pos);

    return 0.5 * vel.z * vel.z + Phi_Rz - Phi_R0;
}

double radial_random_energy(const Vec3& pos, const Vec3& vel)
{
    return random_energy(pos, vel) - vertical_energy(pos, vel);
}

double radial_action(const Vec3& pos, const Vec3& vel)
{
    // Epicyclic approximation:
    // JR ≈ ER / kappa(Rg)
    const double Lz = angular_momentum_z(pos, vel);
    const double Rg = std::abs(Lz) / v_c;

    if (Rg < 1e-12)
        return 0.0;

    const double ER = radial_random_energy(pos, vel);

    if (ER <= 0.0)
        return 0.0;

    const double kappa = std::sqrt(2.0) * v_c / Rg;

    if (kappa <= 0.0)
        return 0.0;

    return ER / kappa;
}

double vertical_action(const Vec3& pos, const Vec3& vel)
{
    // Epicyclic approximation:
    // Jz ≈ Ez / nu, with nu = vc / (q R) near the midplane of the
    // flattened logarithmic potential.

    const double R = std::hypot(pos.x, pos.y);

    if (R < 1e-12 || q <= 0.0 || v_c <= 0.0)
        return 0.0;

    const double Ez = vertical_energy(pos, vel);

    if (Ez <= 0.0)
        return 0.0;

    const double nu = v_c / (q * R);

    if (nu <= 0.0)
        return 0.0;

    return Ez / nu;
}


// ======================================================
// Orbit output
// ======================================================
void simulate_trajectory(Vec3& orbit_pos, Vec3& orbit_vel, const std::string& basename, double& total_time, double& dt, bool use_pert)                   
{
    std::ofstream trajectory_file(basename + "_position.dat");
    std::ofstream parameters_file(basename + "_parameters.dat");
    std::ofstream accel_file     (basename + "_accelerations.dat");
    std::ofstream spiral_file    (basename + "_spiral_profile.dat");

    if (!trajectory_file || !parameters_file || !accel_file || !spiral_file)
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }

    trajectory_file << "# t x y z R phi vx vy vz vR vphi vz_cyl\n";

    parameters_file << "# t Lx Ly Lz "
                    << "KE PE_axis PE_sp PE_total "
                    << "E EJ dE dLz dEJ scale "
                    << "Erand ER Ez JR Jz "
                    << "dErand dER dEz dJR dJz\n";

    accel_file  << "# t ax ay az ax_sp ay_sp az_sp\n";
    spiral_file << "# t R phi z Phi_sp scale Z\n";

    const int steps = static_cast<int>(total_time / dt);
    const int print_interval = std::max(1, steps / 10);

    // Initial total diagnostics.
    const double E0  = total_energy_3D(orbit_pos, orbit_vel, 0.0, use_pert);
    const double Lz0 = angular_momentum_z(orbit_pos, orbit_vel);
    const double EJ0 = jacobi_integral_3D(orbit_pos, orbit_vel, 0.0, use_pert);

    // Initial background heating/action diagnostics.
    // These use the unperturbed logarithmic potential.
    const double Erand0 = random_energy(orbit_pos, orbit_vel);
    const double ER0    = radial_random_energy(orbit_pos, orbit_vel);
    const double Ez0    = vertical_energy(orbit_pos, orbit_vel);
    const double JR0    = radial_action(orbit_pos, orbit_vel);
    const double Jz0    = vertical_action(orbit_pos, orbit_vel);

    for (int i = 0; i <= steps; i++)
    {
        const double sim_time = static_cast<double>(i) * dt;

        Cyl pos_cyl = cart_to_cyl(orbit_pos);
        Cyl vel_cyl = cart_to_cyl_vel(orbit_pos, orbit_vel);

        Vec3 acc_total = total_acceleration_3D(orbit_pos, sim_time, use_pert);
        Vec3 acc_sp    = use_pert ? spiral_acceleration_3D(orbit_pos, sim_time)
                                   : Vec3{0.0, 0.0, 0.0};

        const double KE       = kinetic_energy(orbit_vel);
        const double PE_axis  = potential_energy(orbit_pos);
        const double PE_sp    = use_pert ? spiral_potential_3D(orbit_pos, sim_time) : 0.0;
        const double PE_total = PE_axis + PE_sp;
        const double E        = KE + PE_total;

        Vec3 L = angular_momentum(orbit_pos, orbit_vel);
        const double Lz = L.z;
        const double EJ = E - omega_p * Lz;

        const double dE  = (std::abs(E0)  > 0.0) ? (E  - E0)  / std::abs(E0)  : 0.0;
        const double dLz = (std::abs(Lz0) > 0.0) ? (Lz - Lz0) / std::abs(Lz0) : 0.0;
        const double dEJ = (std::abs(EJ0) > 0.0) ? (EJ - EJ0) / std::abs(EJ0) : 0.0;

        const double scale = use_pert ? perturbation_scaling_factor(sim_time) : 0.0;
        const double Z     = use_pert ? spiral_vertical_envelope(orbit_pos.z)  : 0.0;

        // Background random-energy/action diagnostics.
        // These are not exact integrals during the perturbation. They are
        // osculating background diagnostics used to measure heating.
        const double Erand = random_energy(orbit_pos, orbit_vel);
        const double ER    = radial_random_energy(orbit_pos, orbit_vel);
        const double Ez    = vertical_energy(orbit_pos, orbit_vel);
        const double JR    = radial_action(orbit_pos, orbit_vel);
        const double Jz    = vertical_action(orbit_pos, orbit_vel);

        const double dErand = (std::abs(Erand0) > 0.0) ? (Erand - Erand0) / std::abs(Erand0) : 0.0;
        const double dER    = (std::abs(ER0)    > 0.0) ? (ER    - ER0)    / std::abs(ER0)    : 0.0;
        const double dEz    = (std::abs(Ez0)    > 0.0) ? (Ez    - Ez0)    / std::abs(Ez0)    : 0.0;
        const double dJR    = (std::abs(JR0)    > 0.0) ? (JR    - JR0)    / std::abs(JR0)    : 0.0;
        const double dJz    = (std::abs(Jz0)    > 0.0) ? (Jz    - Jz0)    / std::abs(Jz0)    : 0.0;

        trajectory_file << std::setprecision(10)
                        << sim_time << " "
                        << orbit_pos.x << " " << orbit_pos.y << " " << orbit_pos.z << " "
                        << pos_cyl.R << " " << pos_cyl.phi << " "
                        << orbit_vel.x << " " << orbit_vel.y << " " << orbit_vel.z << " "
                        << vel_cyl.R << " " << vel_cyl.phi << " " << vel_cyl.z << "\n";

        parameters_file << std::setprecision(10)
                        << sim_time << " "
                        << L.x << " " << L.y << " " << L.z << " "
                        << KE << " " << PE_axis << " " << PE_sp << " " << PE_total << " "
                        << E << " " << EJ << " "
                        << dE << " " << dLz << " " << dEJ << " " << scale << " "
                        << Erand << " " << ER << " " << Ez << " " << JR << " " << Jz << " "
                        << dErand << " " << dER << " " << dEz << " " << dJR << " " << dJz << "\n";

        accel_file << std::setprecision(10)
                   << sim_time << " "
                   << acc_total.x << " " << acc_total.y << " " << acc_total.z << " "
                   << acc_sp.x << " " << acc_sp.y << " " << acc_sp.z << "\n";

        spiral_file << std::setprecision(10)
                    << sim_time << " "
                    << pos_cyl.R << " " << pos_cyl.phi << " " << orbit_pos.z << " "
                    << PE_sp << " " << scale << " " << Z << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "[Orbit: " << basename << "] "
                      << "t = " << sim_time
                      << " | E = " << E
                      << " | Lz = " << Lz
                      << " | EJ = " << EJ
                      << " | dE = " << dE
                      << " | dLz = " << dLz
                      << " | dEJ = " << dEJ
                      << " | JR = " << JR
                      << " | Jz = " << Jz
                      << " | dJR = " << dJR
                      << " | dJz = " << dJz
                      << " | scale = " << scale << "\n";
        }

        Leapfrog_integrator_3D(orbit_pos, orbit_vel, dt, sim_time, use_pert);
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
    spiral_file.close();
}

// ======================================================
// Population summary output
// ======================================================
void simulate_population_summary_3D(std::vector<InitialStar3D>& stars, double total_time, double dt, bool use_pert, std::string& output_filename)                                   
{
    std::ofstream out(output_filename);
    out << "# id R0 Rf Rg0 Rgf Lz0 Lzf dLz dRg "
        << "E0 Ef dE EJ0 EJf dEJ "
        << "Erand0 Erandf dErand ER0 ERf dER Ez0 Ezf dEz "
        << "JR0 JRf dJR Jz0 Jzf dJz "
        << "z0 zf vz0 vzf\n";

    const int steps = static_cast<int>(total_time / dt);
    const int progress_interval = std::max(1, static_cast<int>(stars.size()) / 20);
        
    double sum_dLz = 0.0;
    double sum_dRg = 0.0;
    double sum_dJR = 0.0;
    double sum_dJz = 0.0;
    double sum_dLz2 = 0.0;

    for (size_t n = 0; n < stars.size(); n++)
    {
        Vec3 pos = stars[n].pos;
        Vec3 vel = stars[n].vel;

        const double R0  = std::hypot(pos.x, pos.y);
        const double z0  = pos.z;
        const double vz0 = vel.z;

        const double Lz0 = angular_momentum_z(pos, vel);
        const double Rg0 = Lz0 / v_c;

        const double E0  = total_energy_3D(pos, vel, 0.0, use_pert);
        const double EJ0 = jacobi_integral_3D(pos, vel, 0.0, use_pert);

        const double Erand0 = random_energy(pos, vel);
        const double ER0    = radial_random_energy(pos, vel);
        const double Ez0    = vertical_energy(pos, vel);
        const double JR0    = radial_action(pos, vel);
        const double Jz0    = vertical_action(pos, vel);

        for (int i = 0; i < steps; i++)
        {
            const double sim_time = static_cast<double>(i) * dt;
            Leapfrog_integrator_3D(pos, vel, dt, sim_time, use_pert);
        }

        const double Rf  = std::hypot(pos.x, pos.y);
        const double zf  = pos.z;
        const double vzf = vel.z;

        const double Lzf = angular_momentum_z(pos, vel);
        const double Rgf = Lzf / v_c;

        const double Ef  = total_energy_3D(pos, vel, total_time, use_pert);
        const double EJf = jacobi_integral_3D(pos, vel, total_time, use_pert);

        const double Erandf = random_energy(pos, vel);
        const double ERf    = radial_random_energy(pos, vel);
        const double Ezf    = vertical_energy(pos, vel);
        const double JRf    = radial_action(pos, vel);
        const double Jzf    = vertical_action(pos, vel);

        const double dLz    = Lzf - Lz0;
        const double dRg    = Rgf - Rg0;
        const double dE     = Ef - E0;
        const double dEJ    = EJf - EJ0;
        const double dErand = Erandf - Erand0;
        const double dER    = ERf - ER0;
        const double dEz    = Ezf - Ez0;
        const double dJR    = JRf - JR0;
        const double dJz    = Jzf - Jz0;

        out << std::setprecision(10)
            << stars[n].id << " "
            << R0 << " " << Rf << " "
            << Rg0 << " " << Rgf << " "
            << Lz0 << " " << Lzf << " "
            << dLz << " " << dRg << " "
            << E0 << " " << Ef << " " << dE << " "
            << EJ0 << " " << EJf << " " << dEJ << " "
            << Erand0 << " " << Erandf << " " << dErand << " "
            << ER0 << " " << ERf << " " << dER << " "
            << Ez0 << " " << Ezf << " " << dEz << " "
            << JR0 << " " << JRf << " " << dJR << " "
            << Jz0 << " " << Jzf << " " << dJz << " "
            << z0 << " " << zf << " "
            << vz0 << " " << vzf << "\n";

        sum_dLz  += dLz;
        sum_dRg  += dRg;
        sum_dJR  += dJR;
        sum_dJz  += dJz;
        sum_dLz2 += dLz * dLz;

        if (n % progress_interval == 0)
        {
            std::cout << "Integrated " << n << " / " << stars.size() << " stars\n";           
        }
    }

    const double N = static_cast<double>(stars.size());

    if (N > 0.0)
    {
        std::cout << "Mean dLz = " << sum_dLz / N << "\n";
        std::cout << "RMS dLz  = " << std::sqrt(sum_dLz2 / N) << "\n";
        std::cout << "Mean dRg = " << sum_dRg / N << "\n";
        std::cout << "Mean dJR = " << sum_dJR / N << "\n";
        std::cout << "Mean dJz = " << sum_dJz / N << "\n";
    }

    std::cout << "Population summary written to " << output_filename << "\n";           
}

int main(int argc, char** argv)
{
    double total_time = 4000.0; // Myr
    double dt = 0.01;          // Myr
    bool use_pert = true;

    std::string initial_conditions_file = "DF_initial_conditions_3D.dat";
    std::string output_file = "deltaLz_from_DF_3D.dat";

    // -1 means use all stars from the initial-condition file.
    int max_stars = -1;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-pert")
        {
            use_pert = true;
        }
        else if (arg == "-no_pert")
        {
            use_pert = false;
        }
        else if (arg == "-n" && i + 1 < argc)
        {
            max_stars = std::stoi(argv[++i]);
        }
        else if (arg == "-dt" && i + 1 < argc)
        {
            dt = std::stod(argv[++i]);
        }
        else if (arg == "-t" && i + 1 < argc)
        {
            total_time = std::stod(argv[++i]);
        }
        else if (arg == "-ic" && i + 1 < argc)
        {
            initial_conditions_file = argv[++i];
        }
        else if (arg == "-out" && i + 1 < argc)
        {
            output_file = argv[++i];
        }
        else
        {
            std::cerr << "Unknown or incomplete option: " << arg << "\n";
            std::cerr << "Usage: ./LogPotExec [-pert|-no_pert] [-n N] [-dt DT] [-t T] [-ic FILE] [-out FILE]\n";
            return 1;
        }
    }

    try
    {
        std::cout << "Running 3D logarithmic-potential population integration.\n";
        std::cout << "Reading initial conditions from "
                  << initial_conditions_file << "\n";
        std::cout << "Output summary file: " << output_file << "\n";
        std::cout << "total_time = " << total_time
                  << ", dt = " << dt << "\n";
        std::cout << "v_c = " << v_c
                  << ", q = " << q << "\n";
        std::cout << "use_pert = " << std::boolalpha << use_pert
                  << ", m = " << m
                  << ", R_CR = " << R_CR
                  << ", omega_p = " << omega_p
                  << ", pert_strength = " << pert_strength
                  << ", spiral_z_scale = " << spiral_z_scale << "\n";

        std::vector<InitialStar3D> stars = read_initial_conditions_3D(initial_conditions_file, max_stars);
        std::cout << "Loaded " << stars.size() << " stars from initial-conditions file.\n";
                  

        simulate_population_summary_3D(stars, total_time, dt, use_pert, output_file);
        std::cout << "Finished 3D population integration.\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}