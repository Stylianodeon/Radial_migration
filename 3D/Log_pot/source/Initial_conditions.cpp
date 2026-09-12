#include "Initial_conditions.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

InitialConditions3D::InitialConditions3D(int N_R_grid, double R_min_, double R_max_)
                                        : N_grid(N_R_grid), R_min(R_min_), R_max(R_max_)     
{
    precompute_cdf();
}

//Surface density Σ
double InitialConditions3D::Sigma(double R) 
{
    return Sigma_o * std::exp(-R / R_o);
}

//radial velocity dispersion σ_R
double InitialConditions3D::sigma_R(double R) 
{
    return sigma_R0 * std::exp(-(R - R_ref) / R_sigma);
}

//vertical velocity dispersion σ_z
double InitialConditions3D::sigma_z(double R) 
{
    return sigma_z0 * std::exp(-(R - R_ref) / R_sigma);
}

//star frequency at certain radius Ω
double InitialConditions3D::Omega(double R) 
{
    return v_c / R;
}

// Flat rotation curve: kappa = sqrt(Ω)
double InitialConditions3D::kappa(double R)
{
    return std::sqrt(2.0) * Omega(R);
}

double InitialConditions3D::nu(double R) 
{
    // Near-midplane vertical frequency for Φ = 0.5 vc^2 ln(R^2 + z^2/q^2).
    return v_c / (q * R);
}

double InitialConditions3D::E_circ(double R) 
{
    return 0.5 * v_c * v_c + 0.5 * v_c * v_c * std::log(R * R);       
}

double InitialConditions3D::L_circ(double R) 
{
    return R * v_c;
}

//Find the asymmeetric driift to obtain an estimate of the mean azimuthal velocity for the star distribution
double InitialConditions3D::asymmetric_drift(double R) 
{
    const double Om = Omega(R);
    const double kap = kappa(R);

    if (R <= 0.0 || Om <= 0.0 || kap <= 0.0 || v_c <= 0.0)
        return 0.0;

    const double sR = sigma_R(R);
    const double sphi = sR * (kap / (2.0 * Om)); //from epicycle approximation BT08 eqn 8.117

    const double sR2 = sR * sR;
    const double sphi2 = sphi * sphi;

    //calculate asymmetric drift v_α
    const double dlnSigma_dlnR = -R / R_o;
    const double dlnsR2_dlnR = -2.0 * R / R_sigma;
    const double bracket = - dlnSigma_dlnR - dlnsR2_dlnR - (1.0 - sphi2 / sR2);

    double v_a = (0.5 * sR2 / v_c) * bracket;
    return v_a;
}

//draw radii from the exponential disc p(R) ∝ RΣ(R)
void InitialConditions3D::precompute_cdf()
{
    if (N_grid < 2)
    {
        throw std::runtime_error("N_grid must be >= 2");       
    }
    if (R_min <= 0.0 || R_o <= 0.0)
    {
        throw std::runtime_error("InitialConditions3D: R_min and R_o must be positive");     
    }

    R_vals.resize(N_grid);
    CDF_vals.resize(N_grid);

    const double dR = (R_max - R_min) / static_cast<double>(N_grid - 1); //calculate grid spacing
        
    for (int i = 0; i < N_grid; i++)
    {
        R_vals[i] = R_min + static_cast<double>(i) * dR;
            
    }

    CDF_vals[0] = 0.0;
    double total = 0.0;

    for (int i = 1; i < N_grid; i++)
    {
        const double R1 = R_vals[i - 1];
        const double R2 = R_vals[i];

        // dP ∝ 2 π R Sigma(R) dR.
        const double weight1 = R1 * Sigma(R1);
        const double weight2 = R2 * Sigma(R2);

        // Trapezoidal integration of the radial probability density.
        total += 0.5 * (weight1 + weight2) * dR;
        CDF_vals[i] = total;
    }

    if (total <= 0.0 || !std::isfinite(total))
    {
        throw std::runtime_error("InitialConditions3D: invalid CDF total");     
    }

    for (double& value : CDF_vals)
    {
        value /= total;
    }
}

//choose a random R value from the CDF
double InitialConditions3D::sample_R_exponential(std::mt19937& gen)
{
    if (CDF_vals.empty() || R_vals.size() != CDF_vals.size())
    {
        throw std::runtime_error("InitialConditions3D: CDF not initialized correctly");
            
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0); //choose random number between 0 and 1
    const double u = uniform01(gen);

    auto it = std::lower_bound(CDF_vals.begin(), CDF_vals.end(), u);
    int idx = static_cast<int>(std::distance(CDF_vals.begin(), it));

    if (idx < 1)
        idx = 1;

    if (idx >= N_grid)
        idx = N_grid - 1;

    const double R1 = R_vals[idx - 1];
    const double R2 = R_vals[idx];
    const double C1 = CDF_vals[idx - 1];
    const double C2 = CDF_vals[idx];

    const double denominator = C2 - C1;

    if (denominator <= 0.0)
        return R1;

    const double fraction = (u - C1) / denominator;

    return R1 + fraction * (R2 - R1);
}

//Create the initial conditions for one star
SampleOrbit3D InitialConditions3D::sample_star(std::mt19937& gen)
{
    std::uniform_real_distribution<double> uniform_phi(0.0, 2.0 * M_PI); //find a random azimuthal angle
    const int max_attempts = 10000;

    for (int attempt = 0; attempt < max_attempts; attempt++)
    {
        const double R = sample_R_exponential(gen);

        if (!std::isfinite(R) || R <= 0.0)
            continue;

        const double phi = uniform_phi(gen);
        const double Om = Omega(R);
        const double kap = kappa(R);
        const double sR = sigma_R(R);
        const double sz = sigma_z(R);

        if (Om <= 0.0 || kap <= 0.0 || sR <= 0.0 || sz <= 0.0)
            continue;

        const double sphi = sR * (kap / (2.0 * Om));
        const double v_a = asymmetric_drift(R);
        const double mean_vphi = v_c - v_a;
        const double vertical_frequency = nu(R);

        if (vertical_frequency <= 0.0)
            continue;

        // Approximate harmonic-oscillator vertical scale height.
        const double z_scale = sz / vertical_frequency;

        std::normal_distribution<double> normal_vR(0.0, sR);
        std::normal_distribution<double> normal_vphi(mean_vphi, sphi);
        std::normal_distribution<double> normal_vz(0.0, sz);
        std::normal_distribution<double> normal_z(0.0, z_scale);

        const double v_R = normal_vR(gen);
        const double v_phi = normal_vphi(gen);
        const double v_z = normal_vz(gen);
        const double z = normal_z(gen);

        // Reject extreme or invalid initial conditions.
        if (!std::isfinite(v_R) || !std::isfinite(v_phi) || !std::isfinite(v_z) || !std::isfinite(z))  
        {
            continue;
        }

        if (v_phi <= 0.0)
            continue;

        if (std::abs(v_R) > sigma_cut * sR)
            continue;

        if (std::abs(v_z) > sigma_cut * sz)
            continue;

        if (std::abs(z) > max_abs_z)
            continue;

        Cyl pos_cyl{R, phi, z};
        Cyl vel_cyl{v_R, v_phi, v_z};

        Vec3 pos = cyl_to_cart(pos_cyl);
        Vec3 vel = cyl_to_cart_vel(pos_cyl, vel_cyl);

        const double Lz = R * v_phi;
        const double Rg = Lz / v_c;

        return {R, phi, z, v_R, v_phi, v_z, Lz, Rg, pos, vel};    
    }
    throw std::runtime_error("InitialConditions3D: failed to sample a valid star after many attempts");
}

void InitialConditions3D::write_initial_conditions(const std::string& filename, int N_samples, unsigned int seed)
{
    if (N_samples <= 0)
    {
        throw std::runtime_error("InitialConditions3D: N_samples must be positive");        
    }

    std::mt19937 gen(seed);
    std::ofstream out(filename);

    if (!out)
    {
        throw std::runtime_error("InitialConditions3D: could not open output file " + filename);        
    }

    out << std::setprecision(10);
    out << "# id R phi z v_R v_phi v_z Lz Rg "
        << "x y z_cart vx vy vz_cart\n";

    for (int i = 0; i < N_samples; i++)
    {
        const SampleOrbit3D star = sample_star(gen);

        out << i << " " << star.R << " " << star.phi << " " << star.z << " " << star.v_R << " " << star.v_phi << " "
            << star.v_z << " " << star.Lz << " " << star.Rg << " " << star.pos.x << " " << star.pos.y << " " << star.pos.z << " "
            << star.vel.x << " " << star.vel.y << " " << star.vel.z << "\n";
    }

    out.close();
    std::cout << "Generated " << N_samples << " 3D warm-disc initial conditions in " << filename << "\n";             
}

/*
int main()
{
   
    int N_samples = 30000;
    unsigned int seed = 2;

    InitialConditions3D ic(
        30000,  // radial CDF grid resolution
        3.0,    // minimum radius[kpc]
        15.0);  // maximum  radius [kpc]

    ic.write_initial_conditions("DF_initial_conditions_3D.dat", N_samples, seed);   
    
    return 0;
}
*/
