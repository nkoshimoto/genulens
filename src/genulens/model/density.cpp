#include "genulens/model/density.hpp"

#include <algorithm>
#include <cmath>

namespace genulens::model {

ComponentDensities approximate_density(double distance_pc, GalacticCoordinates coordinates)
{
    ComponentDensities densities;
    const double z = distance_pc * std::sin(coordinates.b_deg * 3.1415926535897932385 / 180.0);
    const double disk = std::exp(-std::abs(z) / 300.0) * std::exp(-distance_pc / 12000.0);
    for (int i = 0; i < 7; ++i) {
        densities.values[static_cast<std::size_t>(i)] = disk * (0.01 + 0.004 * i);
    }
    densities.values[7] = 0.003 * std::exp(-std::abs(z) / 900.0);
    densities.values[8] = 0.15 * std::exp(-std::pow((distance_pc - 8000.0) / 2500.0, 2.0));
    densities.values[9] = 0.0;
    densities.values[10] = 1e-4 * std::exp(-distance_pc / 20000.0);
    return densities;
}

} // namespace genulens::model

