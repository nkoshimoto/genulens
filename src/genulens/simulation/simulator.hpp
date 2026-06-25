#pragma once

#include "genulens/model/galactic_model.hpp"
#include "genulens/rng.hpp"
#include "genulens/simulation/likelihood.hpp"

namespace genulens {

class EventSimulator {
public:
    explicit EventSimulator(GenulensConfig config);

    SimulationResult simulate(LikelihoodFunction likelihood = {});

private:
    GenulensConfig config_;
};

SimulationResult simulate(const GenulensConfig &config, LikelihoodFunction likelihood = {});

} // namespace genulens

