#include "Dehnen_DF.h"
#include "VecUtils.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>


double WarmDiskDF::Sigma(double R) const 
{
    return Sigma_o * std::exp(-R / R_o);
}

void WarmDiskDF::precompute_cdf()
{
    if (N_grid < 2 || !std::isfinite(R_min) || !std::isfinite(R_max) || R_min <= 0.0 || R_max <= R_min)
        throw std::runtime_error("WarmDiskDF: require N_grid >= 2 and 0 < Rg_min < Rg_max");

    R_vals.resize(N_grid);
    CDF_vals.resize(N_grid);
    const double dR = (R_max - R_min) / (N_grid - 1);
    double total = 0.0;
    double previous_weight = 0.0;
    for (int i = 0; i < N_grid; i++)
    {
        const double Rg = R_min + i * dR;
        R_vals[i] = Rg;
        const double sigma = sigma_R(Rg);
        const long double beta = (static_cast<long double>(v_c) / sigma) * (v_c / sigma);

        if (!(sigma > 0.0) || !std::isfinite(beta) || beta <= 1.0L)
        {
            throw std::runtime_error("WarmDiskDF: the logarithmic Shu sampler requires 0 < sigma_R(Rg) < vc");
        }
        // Shu f(E,Lz) = [2 Omega/kappa] Sigma(Rg)/(2 pi sigma^2) * exp[-(E-Ec(Lz))/sigma^2], for prograde Lz.           
        const long double shape = (beta - 1.0L) / 2.0L;
        const long double log_I = beta / 2.0L - std::log(2.0L) + std::lgamma(shape) - shape * std::log(beta / 2.0L);                       
        const double weight = static_cast<double>(Rg * Sigma(Rg) * std::exp(log_I) / sigma);

        if (!std::isfinite(weight) || weight <= 0.0)
        {
            throw std::runtime_error("WarmDiskDF: invalid guiding-radius weight");
        }

        if (i > 0) 
        {
            total += 0.5 * (previous_weight + weight) * dR;
        }
        CDF_vals[i] = total;
        previous_weight = weight;
    }
    if (!std::isfinite(total) || total <= 0.0)
    {
        throw std::runtime_error("WarmDiskDF: invalid CDF normalization");
    }
    for (double& value : CDF_vals) 
    {
        value /= total;
    }
}

double WarmDiskDF::sample_R(std::mt19937& gen) 
{
    if (CDF_vals.empty() || R_vals.size() != CDF_vals.size()) 
    {
        throw std::runtime_error("Distribution values not initialized correclty");
    }
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double u = dist(gen);

    // Binary search for u in CDF
    auto it = std::lower_bound(CDF_vals.begin(), CDF_vals.end(), u);
    int idx = static_cast<int>(std::distance(CDF_vals.begin(), it));
    if (idx < 1) idx = 1;
    if (idx >= N_grid) idx = N_grid - 1;

    double R1 = R_vals[idx - 1];
    double R2 = R_vals[idx];
    double C1 = CDF_vals[idx - 1];
    double C2 = CDF_vals[idx];
    const double denom = (C2 - C1);
    if (denom <= 0.0) 
    {
        // fallback: return left edge of the interval
        return R1;
    }
    double t = (u - C1) / denom;
    return R1 + t * (R2 - R1);
}

double WarmDiskDF::E_circ(double R) const 
{
    return 0.5 * v_c * v_c + 0.5 * v_c * v_c * std::log(R * R);
}

double WarmDiskDF::L_circ(double R) const 
{
    return R * v_c;
}

double WarmDiskDF::Omega(double R) const 
{
    return v_c / R;
}

// Constructor to initialize values
WarmDiskDF::WarmDiskDF(int N_grid_, double R_min_, double R_max_) : N_grid(N_grid_), R_min(R_min_), R_max(R_max_)
{
    precompute_cdf();
}

double WarmDiskDF::sigma_R(double R) const
{
    return sigma_o * std::exp(-R / R_sigma);
}

SampleOrbit WarmDiskDF::sample_orbit(std::mt19937& gen)
{
    const double Rg = sample_R(gen);
    const double sigma = sigma_R(Rg);
    const double beta = (v_c / sigma) * (v_c / sigma);

    // At fixed Lz=vc*Rg, vR is Gaussian and independent of R. For the
    // logarithmic potential, u=(Rg/R)^2 has an exact gamma distribution:
    // p(u|Rg) proportional to u^((beta-3)/2) exp(-beta*u/2).
    std::gamma_distribution<double> radial_ratio((beta - 1.0) / 2.0, 2.0 / beta);
    std::normal_distribution<double> radial_velocity(0.0, sigma);
    const double u = radial_ratio(gen);
    if (!std::isfinite(u) || u <= 0.0)
    {
        throw std::runtime_error("WarmDiskDF: invalid gamma sample");
    }
    const double R = Rg / std::sqrt(u);
    const double L = L_circ(Rg);
    const double v_R = radial_velocity(gen);
    const double v_phi = L / R;
    return {R, L, v_R, v_phi};
}

double WarmDiskDF::kappa(double R) const
{
    return std::sqrt(2.0) * Omega(R);
}

#ifdef BUILD_SHU_DF_MAIN
int main()
{
    std::mt19937 gen(42); // Reproducible initial conditions.
    WarmDiskDF df;

    const int N_samples = 30000;

    // Open .dat file (space-separated)
    std::ofstream DF("DF_initial_conditions.dat");
    if (!DF) 
    {
        std::cerr << "Error: could not open DF_initial_conditions.dat for writing.\n";
        return 1;
    }

    DF << std::setprecision(17);
    DF << "# R  L  v_R  v_phi   x  y\n";
    std::uniform_real_distribution<double> uniform_phi(0.0, 2.0 * M_PI); //generate random φ values for a specific R coordinate

    for (int i = 0; i < N_samples; i++)
    {
        SampleOrbit orb = df.sample_orbit(gen);
        double phi = uniform_phi(gen);
        Cyl cyl_coords{orb.R, phi};
        Vec2 pos = cyl_to_cart(cyl_coords);

        DF << orb.R << " " << orb.L << " " << orb.v_R << " " << orb.v_phi <<  " "  << pos.x << " " << pos.y << "\n"; 
    }

    DF.close();
    std::cout << "Generated " << N_samples << " orbits based on Shu Distribution Function\n";
    return 0;
}
#endif


