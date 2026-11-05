#include "SpiralPerturbation.h"
#include "LogPot.h"

#include <cmath>

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

void spiral_amplitude(double R, double scale, double& A, double& dA)
{
    const double PhiS = (2.0 * M_PI * G * e_s * Sigma_o * std::exp(-R / R_o)) / (alpha / R);
    const double W = spiral_radial_envelope(R);

    A = pert_strength * scale * PhiS * W;
    dA = pert_strength * scale * PhiS * ((1.0 / R - 1.0 / R_o) * W + d_spiral_radial_envelope_dR(R));

}

double spiral_phase(double R, double phi, double time)
{
    return alpha * std::log(R / R_CR) + m * omega_p * time - m * phi;
}

double spiral_potential(const Cyl& pos, double time, double scale)
{
    double A, dA;
    spiral_amplitude(pos.R, scale, A, dA);
    return A * std::cos(spiral_phase(pos.R, pos.phi, time));
}

Vec2 spiral_gradient(const Vec2& pos, double time, double scale)
{
    const double R = std::hypot(pos.x, pos.y);
    double A, dA;
    spiral_amplitude(R, scale, A, dA);
    const double phase = spiral_phase(R, std::atan2(pos.y, pos.x), time);
    const double radial = dA * std::cos(phase) - A * std::sin(phase) * alpha / R;
    const double azimuthal = m * A * std::sin(phase) / R;
    return {(radial * pos.x - azimuthal * pos.y) / R, (radial * pos.y + azimuthal * pos.x) / R};

}
