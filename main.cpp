#include "concave_hull.hpp"
#include "matplotlibcpp.hpp"

#include <chrono>
#include <random>
#include <vector>

namespace plt = matplotlibcpp;

int main()
{
    constexpr int number_of_points = 1000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-100.0, 100.0);

    std::vector<double> coords;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    for (int i = 0; i < number_of_points; ++i)
    {
        const double x = dist(gen);
        const double y = dist(gen);
        coords.push_back(x);
        coords.push_back(y);
        x_coords.push_back(x);
        y_coords.push_back(y);
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    const auto hull_indices = geometry::constructConcaveHull(coords, 0.3);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time (s): " << (t2 - t1).count() / 1e9 << std::endl;

    std::vector<double> hull_coords = delaunator::getHullCoordinates(hull_indices, coords);

    std::vector<double> x_hull;
    std::vector<double> y_hull;

    for (std::size_t i = 0; i < hull_coords.size(); i += 2)
    {
        const double x = hull_coords[i];
        const double y = hull_coords[i + 1];

        x_hull.push_back(x);
        y_hull.push_back(y);
    }

    plt::plot(x_hull, y_hull, {{"marker", "x"}, {"markersize", "3"}});

    plt::plot(x_coords, y_coords, {{"marker", "."}, {"markersize", "1"}, {"linestyle", "None"}});

    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}