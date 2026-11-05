#include "EffectivePotential.h"
#include "LogPot.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

using std::cos;
using std::exp;
using std::log;
using std::string;


static constexpr double minimum_radius = 1e-6; // The logarithmic model is singular at R=0. Exclude this central region (kpc)


static void spiral_amplitude(double R, double scale, double& A, double& dA)
{
    const double PhiS = (2.0 * M_PI * G * e_s * Sigma_o * std::exp(-R / R_o)) / (alpha / R);
    const double W = spiral_radial_envelope(R);

    A = pert_strength * scale * PhiS * W;
    dA = pert_strength * scale * PhiS * ((1.0 / R - 1.0 / R_o) * W + d_spiral_radial_envelope_dR(R));

}

static Vec2 effective_gradient(const Vec2& pos, bool use_pert, double time, double scale)
{
    const double R = std::hypot(pos.x, pos.y);

    double dR = v_c * v_c / R - omega_p * omega_p * R;
    double dphi = 0.0;
    if (use_pert && pert_strength != 0.0 && scale != 0.0)
    {
        double A, dA;
        spiral_amplitude(R, scale, A, dA);
        const double phase = alpha * std::log(R / R_CR) + m * omega_p * time - m * std::atan2(pos.y, pos.x);
        dR += dA * std::cos(phase) - A * std::sin(phase) * alpha / R;
        dphi = m * A * std::sin(phase) / R;
    }
    return {(dR * pos.x - dphi * pos.y) / R, (dR * pos.y + dphi * pos.x) / R};
}

static double radial_equilibrium(double R, double sign, double scale)
{
    double A, dA;
    spiral_amplitude(R, scale, A, dA);
    return v_c * v_c / R - omega_p * omega_p * R + sign * dA;
}

static double bisect_equilibrium(double lo, double hi, double sign, double scale)
{
    double flo = radial_equilibrium(lo, sign, scale);
    for (int i = 0; i < 100; i++)
    {
        const double mid = lo + 0.5 * (hi - lo);
        const double fm = radial_equilibrium(mid, sign, scale);
        if (!std::isfinite(fm))
            throw std::runtime_error("Non-finite radial equilibrium equation");
        if (fm == 0.0 || hi - lo <= 1e-12 * std::max(1.0, mid))
            return mid;
        if (std::signbit(flo) != std::signbit(fm))
            hi = mid;
        else
        {
            lo = mid;
            flo = fm;
        }
    }
    throw std::runtime_error("Lagrange-point bisection did not converge");
}

double effective_potential(const Vec2 &cart_pos, const bool &use_pert)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);
    if (cyl_pos.R <= minimum_radius)
        return std::numeric_limits<double>::quiet_NaN();

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

Vec2 effective_potential_gradient(const Vec2& pos, const bool& use_pert)
{
    return effective_gradient(pos, use_pert, 0.0, 1.0);
}

void scan_Lagrange_points(double xmin, double xmax, double ymin, double ymax, double step, const string& filename, bool use_pert, double sim_time)
{
    const double nx = std::floor((xmax - xmin) / step + 1e-10);
    const double ny = std::floor((ymax - ymin) / step + 1e-10);
    std::ofstream fout(filename);
    std::ofstream points("lagrange_points.dat");
    if (!fout || !points)
        throw std::runtime_error("Cannot open Lagrange-point output files");
    fout << std::setprecision(17);
    points << std::setprecision(17);
    fout << "# x y |grad(Phi_eff)|; inertial coordinates; snapshot_time_Myr=" << sim_time << "\n";
    points << "# Instantaneous stationary points; inertial coordinates; snapshot_time_Myr=" << sim_time
           << "\n# x y |grad(Phi_eff)|\n";
    const double scale = use_pert ? perturbation_scaling_factor(sim_time) : 0.0;

    // The grid is for plotting only; its spacing does not control the root solver.
    for (std::size_t j = 0; j <= static_cast<std::size_t>(ny); ++j)
    {
        const double y = ymin + j * step;
        for (std::size_t i = 0; i <= static_cast<std::size_t>(nx); ++i)
        {
            const double x = xmin + i * step;
            const Vec2 g = effective_gradient({x, y}, use_pert, sim_time, scale);
            fout << x << ' ' << y << ' ' << std::hypot(g.x, g.y) << '\n';
        }
        fout << '\n';
    }

    const double rmin = std::max(minimum_radius, std::hypot(std::clamp(0.0, xmin, xmax), std::clamp(0.0, ymin, ymax)));
    const double rmax = std::hypot(std::max(std::abs(xmin), std::abs(xmax)), std::max(std::abs(ymin), std::abs(ymax)));

    if (scale == 0.0 || pert_strength == 0.0 || e_s == 0.0 || Sigma_o == 0.0 || G == 0.0)
    {
        if (omega_p != 0.0)
        {
            const double radius = v_c / std::abs(omega_p);
            points << "# Axisymmetric equilibrium ring: R=" << radius << " kpc; no isolated Lagrange points.\n";

            if (radius < rmin || radius > rmax) points << "# The equilibrium ring lies outside the search bounds.\n";
        }
        else
            points << "# No finite-radius equilibrium: pattern speed is zero.\n";
        return;
    }
    if (rmax <= rmin)
        return;

    // Resolve sign-changing roots of both radial branches on a logarithmic mesh.
    // This resolution is independent of the Cartesian plotting grid. Degenerate
    // double roots at a bifurcation require a separate continuation analysis.
    constexpr int radial_intervals = 16384;
    const double log_span = std::log(rmax / rmin);
    std::size_t count = 0;
    for (int parity = 0; parity < 2; ++parity)
    {
        const double sign = parity == 0 ? 1.0 : -1.0;
        double lo = rmin;
        double flo = radial_equilibrium(lo, sign, scale);
        double previous_root = -1.0;
        for (int i = 1; i <= radial_intervals; ++i)
        {
            const double hi = i == radial_intervals ? rmax : rmin * std::exp(log_span * i / radial_intervals);
            const double fhi = radial_equilibrium(hi, sign, scale);
            if (!std::isfinite(flo) || !std::isfinite(fhi))
                throw std::runtime_error("Non-finite radial equilibrium equation");
            double root = -1.0;
            if (flo == 0.0) root = lo;
            else if (fhi == 0.0) root = hi;
            else if (std::signbit(flo) != std::signbit(fhi))
                root = bisect_equilibrium(lo, hi, sign, scale);
            if (root > minimum_radius &&
                std::abs(root - previous_root) > 1e-10 * std::max(1.0, root))
            {
                previous_root = root;
                for (int arm = 0; arm < static_cast<int>(m); ++arm)
                {
                    const int n = 2 * arm + parity;
                    const double phi = (alpha * std::log(root / R_CR) - n * M_PI) / m + omega_p * sim_time;

                    const Vec2 pos{root * std::cos(phi), root * std::sin(phi)};
                    if (pos.x < xmin || pos.x > xmax || pos.y < ymin || pos.y > ymax)
                        continue;
                    const Vec2 g = effective_gradient(pos, true, sim_time, scale);
                    const double residual = std::hypot(g.x, g.y);
                    const double force_scale = v_c * v_c / root + omega_p * omega_p * root;
                    if (!std::isfinite(residual) || residual > 1e-10 * force_scale)
                        throw std::runtime_error("Lagrange-point residual exceeds tolerance");
                    points << pos.x << ' ' << pos.y << ' ' << residual << '\n';
                    ++count;
                }
            }
            lo = hi;
            flo = fhi;
        }
    }
    if (count == 0)
        points << "# No isolated stationary points found within the search bounds.\n";
}

double effective_potential_time(const Vec2 &cart_pos, const bool &use_pert, double sim_time)
{
    Cyl cyl_pos = cart_to_cyl(cart_pos);
    if (cyl_pos.R <= minimum_radius)
        return std::numeric_limits<double>::quiet_NaN();

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

Vec2 effective_potential_gradient_time(const Vec2& pos, const bool& use_pert, double sim_time, double /*total_time*/)
{
    return effective_gradient(pos, use_pert, sim_time, perturbation_scaling_factor(sim_time));
}
