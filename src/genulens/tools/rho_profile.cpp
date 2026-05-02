#include "genulens/tools/rho_profile.hpp"

namespace genulens::tools {

std::vector<RhoProfileRow> rho_profile(const GalacticModel &model, double d_min, double d_max, double d_step)
{
    std::vector<RhoProfileRow> rows;
    if (d_step <= 0.0) return rows;
    for (double distance = d_min; distance <= d_max; distance += d_step) {
        rows.push_back({distance, model.density_at(distance)});
    }
    return rows;
}

} // namespace genulens::tools

