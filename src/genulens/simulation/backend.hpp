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
    std::string run_cli_capture(const std::vector<std::string> &args) const;

private:
    std::vector<std::string> build_event_args(const GenulensConfig &config) const;
    SimulationResult parse_verbosity3_events(const std::string &output, LikelihoodFunction likelihood) const;
};

} // namespace genulens
