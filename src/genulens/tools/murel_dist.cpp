#include "genulens/tools/murel_dist.hpp"

#include <cmath>

namespace genulens::tools {

std::vector<MurelDistributionRow> murel_distribution(RandomEngine &rng, int n_simu, double mu_max, double d_mu)
{
    std::vector<MurelDistributionRow> rows;
    if (n_simu <= 0 || mu_max <= 0.0 || d_mu <= 0.0) return rows;
    const int bins = static_cast<int>(mu_max / d_mu);
    rows.assign(static_cast<std::size_t>(bins), {});
    for (int i = 0; i < bins; ++i) rows[static_cast<std::size_t>(i)].murel = (i + 0.5) * d_mu;
    for (int i = 0; i < n_simu; ++i) {
        const double mu = std::abs(8.0 + 3.0 * rng.gaussian());
        const int bin = static_cast<int>(mu / d_mu);
        if (bin >= 0 && bin < bins) rows[static_cast<std::size_t>(bin)].probability += 1.0;
    }
    for (auto &row : rows) row.probability /= (n_simu * d_mu);
    return rows;
}

} // namespace genulens::tools

