#include "EffectivePotential.h"
#include "LogPot.h"
#include "SpiralPerturbation.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Diagnostic snapshots only: this does not integrate stellar orbits.
int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: LagrangeSnapshotsExec NEW_OUTPUT_DIRECTORY\n";
        return 1;
    }
    const auto original = std::filesystem::current_path();
    try
    {
        const auto output = std::filesystem::absolute(argv[1]);
        if (std::filesystem::exists(output))
            throw std::runtime_error("Choose a new output directory to preserve existing results");
        std::filesystem::create_directories(output);
        std::filesystem::current_path(output);

        const double start = 1900.0, end = 3100.0, cadence = 10.0; // Myr
        const double extent = 12.0, spacing = 0.15;             // kpc
        std::ofstream metadata("metadata.dat"), frames("frames.dat"), spiral("spiral.dat");
        if (!metadata || !frames || !spiral)
            throw std::runtime_error("Cannot open snapshot data files");
        metadata << std::setprecision(17) << "# omega_p v_c extent spacing\n"
                 << omega_p << ' ' << v_c << ' ' << extent << ' ' << spacing << '\n';
        frames << std::setprecision(17) << "# frame time_Myr envelope\n";
        spiral << std::setprecision(17) << "# x y Phi_spiral_full_strength_at_phase_zero\n";
        const int cells = static_cast<int>(std::lround(2.0 * extent / spacing));
        for (int j = 0; j <= cells; ++j)
        {
            for (int i = 0; i <= cells; ++i)
            {
                const Vec2 p{-extent + i * spacing, -extent + j * spacing};
                spiral << p.x << ' ' << p.y << ' '
                       << effective_potential(p, true) - effective_potential(p, false) << '\n';
            }
        }
        const int last = static_cast<int>(std::lround((end - start) / cadence));
        for (int frame = 0; frame <= last; ++frame)
        {
            const double time = start + frame * cadence;
            // The plotting grid is unused here; the radial solver is independent of it.
            scan_Lagrange_points(-extent, extent, -extent, extent, extent,
                                 "scratch_gradient.dat", true, time);
            std::ostringstream name;
            name << "points_" << std::setw(4) << std::setfill('0') << frame << ".dat";
            std::filesystem::rename("lagrange_points.dat", name.str());
            frames << frame << ' ' << time << ' ' << perturbation_scaling_factor(time) << '\n';
        }
        std::filesystem::remove("scratch_gradient.dat");
        metadata.close(); frames.close(); spiral.close();
        if (!metadata || !frames || !spiral)
            throw std::runtime_error("Failed to write snapshot data");
        std::filesystem::current_path(original);
        std::cout << "Wrote " << last + 1 << " equilibrium snapshots to " << output << '\n';
    }
    catch (const std::exception& error)
    {
        std::filesystem::current_path(original);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
