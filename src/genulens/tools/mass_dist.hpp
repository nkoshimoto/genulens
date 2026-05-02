#pragma once

#include "genulens/model/galactic_model.hpp"

#include <vector>

namespace genulens::tools {

struct MassDistributionRow {
    double log_mass = 0.0;
    double dndlogm = 0.0;
};

std::vector<MassDistributionRow> mass_distribution(const GalacticModel &model, int n_mass);

} // namespace genulens::tools

