#pragma once

#include "genulens/model/galactic_model.hpp"

#include <vector>

namespace genulens::tools {

struct RhoProfileRow {
    double distance_pc = 0.0;
    ComponentDensities density;
};

std::vector<RhoProfileRow> rho_profile(const GalacticModel &model, double d_min, double d_max, double d_step);

} // namespace genulens::tools

