#include "genulens/tools/mass_dist.hpp"

#include <cmath>

namespace genulens::tools {

std::vector<MassDistributionRow> mass_distribution(const GalacticModel &model, int n_mass)
{
    std::vector<MassDistributionRow> rows;
    if (n_mass <= 0) return rows;
    const double lo = -3.0;
    const double hi = std::log10(120.0);
    const double step = (hi - lo) / n_mass;
    rows.reserve(static_cast<std::size_t>(n_mass + 1));
    for (int i = 0; i <= n_mass; ++i) {
        const double log_mass = lo + step * i;
        rows.push_back({log_mass, model.imf(std::pow(10.0, log_mass))});
    }
    return rows;
}

} // namespace genulens::tools

