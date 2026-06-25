#pragma once

#include "genulens/options.hpp"
#include "genulens/simulation/likelihood.hpp"
#include "genulens/types.hpp"

#include <string>
#include <vector>

namespace genulens {

class SimulationBackend {
public:
    SimulationResult simulate(const GenulensConfig &config, LikelihoodFunction likelihood = {}) const;
};

} // namespace genulens
