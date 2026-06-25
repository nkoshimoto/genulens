#pragma once

#include "genulens/options.hpp"
#include "genulens/simulation/likelihood.hpp"
#include "genulens/simulation/run_context.hpp"
#include "genulens/types.hpp"

namespace genulens {

class Sampler {
public:
    int run_cli(RunContext &context, int argc, char **argv);
    SimulationResult simulate(RunContext &context, int argc, char **argv,
                              LikelihoodFunction likelihood = {});
    SimulationResult simulate(RunContext &context, const GenulensConfig &config,
                              int argc, char **argv,
                              LikelihoodFunction likelihood = {});
};

} // namespace genulens
