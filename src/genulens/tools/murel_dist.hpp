#pragma once

#include "genulens/model/galactic_model.hpp"
#include "genulens/rng.hpp"

#include <vector>

namespace genulens::tools {

struct MurelDistributionRow {
    double murel = 0.0;
    double probability = 0.0;
};

std::vector<MurelDistributionRow> murel_distribution(RandomEngine &rng, int n_simu, double mu_max, double d_mu);

} // namespace genulens::tools

