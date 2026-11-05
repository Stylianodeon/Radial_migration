#include "EffectivePotential.h"
#include "LogPot.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

static std::vector<Vec2> read_points()
{
    std::ifstream file("lagrange_points.dat");
    require(bool(file), "Missing Lagrange-point output");
    std::vector<Vec2> points;
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line.front() == '#') continue;
        Vec2 p;
        double residual;
        std::istringstream row(line);
        require(bool(row >> p.x >> p.y >> residual), "Malformed point row");
        require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(residual), "Non-finite point");
        require(residual < 1e-12, "Root residual is too large");
        points.push_back(p);
    }
    return points;
}

static std::string point_file()
{
    std::ifstream file("lagrange_points.dat");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

static void check_gradients()
{
    const double h = 1e-5;
    for (Vec2 p : {Vec2{2.0, 1.0}, Vec2{-5.0, 3.0}, Vec2{8.0, -6.0}})
    {
        for (bool pert : {false, true})
        {
            const Vec2 g = effective_potential_gradient(p, pert);
            const double gx = (effective_potential({p.x+h,p.y},pert)-effective_potential({p.x-h,p.y},pert))/(2*h);
            const double gy = (effective_potential({p.x,p.y+h},pert)-effective_potential({p.x,p.y-h},pert))/(2*h);
            require(std::hypot(g.x-gx,g.y-gy)<1e-10, "Static analytic gradient disagrees with potential");
            for (double t : {0.0, 2250.0, 2500.0, 2750.0, 3000.0})
            {
                const Vec2 gt = effective_potential_gradient_time(p, pert, t, 5000.0);
                const double tx = (effective_potential_time({p.x+h,p.y},pert,t)-effective_potential_time({p.x-h,p.y},pert,t))/(2*h);
                const double ty = (effective_potential_time({p.x,p.y+h},pert,t)-effective_potential_time({p.x,p.y-h},pert,t))/(2*h);
                require(std::hypot(gt.x-tx,gt.y-ty)<1e-10, "Timed analytic gradient disagrees with potential");
                const Vec2 force = pert ? total_acceleration(p,t) : LogPot_acc(p);
                require(std::hypot(gt.x+force.x+omega_p*omega_p*p.x,
                                   gt.y+force.y+omega_p*omega_p*p.y)<1e-14,
                        "Effective gradient disagrees with orbit force and centrifugal term");
            }
        }
    }
    require(std::isnan(effective_potential({0.0,0.0},false)), "Origin potential must be masked");
    require(std::isnan(effective_potential_gradient({0.0,0.0},false).x), "Origin must not be an equilibrium");
}

static void check_finder()
{
    std::vector<Vec2> reference;
    for (double step : {0.2, 0.1, 0.05})
    {
        scan_Lagrange_points(-20,20,-20,20,step,"gradient.dat",true,2500.0);
        const auto points = read_points();
        require(points.size()==8, "Expected eight isolated equilibria at peak amplitude");
        for (std::size_t i=0; i<points.size(); ++i)
        {
            const auto p = points[i];
            const double expected = i<4 ? 9.998741834242477 : 10.00125758207147;
            require(std::abs(std::hypot(p.x,p.y)-expected)<1e-9, "Unexpected equilibrium radius");
            if (!reference.empty())
                require(std::hypot(p.x-reference[i].x,p.y-reference[i].y)<1e-12, "Roots depend on plotting grid");
            // Independently check the actual snapshot potential at each root.
            const double h=1e-4;
            const double gx=(effective_potential_time({p.x+h,p.y},true,2500)-effective_potential_time({p.x-h,p.y},true,2500))/(2*h);
            const double gy=(effective_potential_time({p.x,p.y+h},true,2500)-effective_potential_time({p.x,p.y-h},true,2500))/(2*h);
            require(std::hypot(gx,gy)<1e-10, "Reported point is not stationary in snapshot potential");
        }
        reference=points;
    }
    // Equal envelope amplitudes on either side of the peak; only orientation changes.
    scan_Lagrange_points(-20,20,-20,20,1,"gradient.dat",true,2250);
    const auto earlier=read_points();
    scan_Lagrange_points(-20,20,-20,20,1,"gradient.dat",true,2750);
    const auto later=read_points();
    require(earlier.size()==8 && later.size()==8, "Expected eight roots during the ramp");
    const double angle=omega_p*500;
    for (std::size_t i=0; i<earlier.size(); ++i)
    {
        const Vec2 p=earlier[i];
        const Vec2 rotated{p.x*std::cos(angle)-p.y*std::sin(angle),p.x*std::sin(angle)+p.y*std::cos(angle)};
        require(std::hypot(rotated.x-later[i].x,rotated.y-later[i].y)<1e-10, "Incorrect snapshot rotation");
    }
    scan_Lagrange_points(0,20,0,20,1,"gradient.dat",true,2500);
    const auto quadrant=read_points();
    const auto expected=std::count_if(reference.begin(),reference.end(),[](Vec2 p){return p.x>=0 && p.y>=0;});
    require(quadrant.size()==static_cast<std::size_t>(expected), "Incorrect clipping to search bounds");
    for (const auto p:quadrant) require(p.x>=0 && p.y>=0, "Point outside search bounds");
    for (bool pert : {false,true})
    {
        scan_Lagrange_points(-20,20,-20,20,1,"gradient.dat",pert,0);
        require(read_points().empty(), "Axisymmetric ring must not become isolated points");
        require(point_file().find("equilibrium ring: R=10")!=std::string::npos, "Missing corotation ring metadata");
    }
    const double strength=pert_strength;
    pert_strength=0;
    scan_Lagrange_points(-20,20,-20,20,1,"gradient.dat",true,2500);
    require(read_points().empty(), "Zero-strength spiral must produce a ring");
    pert_strength=strength;
    const double speed=omega_p;
    omega_p=0;
    scan_Lagrange_points(-20,20,-20,20,1,"gradient.dat",false,2500);
    require(read_points().empty() && point_file().find("No finite-radius equilibrium")!=std::string::npos,
            "Zero-speed axisymmetric model must have no ring");
    omega_p=speed;
}

int main()
{
    const auto original=std::filesystem::current_path();
    const auto directory=std::filesystem::temp_directory_path()/
        ("lagrange-tests-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try
    {
        std::filesystem::create_directory(directory);
        std::filesystem::current_path(directory);
        check_gradients();
        check_finder();
        std::filesystem::current_path(original);
        std::filesystem::remove_all(directory);
        std::cout << "Effective-potential and Lagrange-point checks passed\n";
    }
    catch(const std::exception& error)
    {
        std::filesystem::current_path(original);
        std::cerr << error.what() << "\nTest artifacts: " << directory << '\n';
        return 1;
    }
}
