#include "OrbitDiagnostics.h"
#include "LogPot.h"

#include <cmath>

using std::log;

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
