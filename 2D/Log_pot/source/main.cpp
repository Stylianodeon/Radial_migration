#include "LogPot.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::cos;
using std::log;

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
            std::ofstream out("sim_results.dat");
            if (!out)
            {
                std::cerr << "Could not open sim_results.dat for writing\n";
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
            std::cout << "Wrote sim_results.dat (6-col DF input)\n";
        }
    }

    {
        std::ifstream fin2("sim_results.dat");
        if (!fin2)
        {
            std::cerr << "Could not open sim_results.dat for RMS calc\n";
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
