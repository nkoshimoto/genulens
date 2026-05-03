#include "genulens/simulation/simulator.hpp"

#include "genulens/simulation/backend.hpp"

#include <algorithm>
#include <cmath>

namespace genulens {

EventSimulator::EventSimulator(GenulensConfig config)
    : config_(std::move(config))
{
}

SimulationResult EventSimulator::simulate(LikelihoodFunction likelihood)
{
    if (!likelihood) {
        likelihood = default_likelihood(config_.observed_tE, config_.observed_tE_error);
    }
    return SimulationBackend().simulate(config_, std::move(likelihood));
}

SimulationResult simulate(const GenulensConfig &config, LikelihoodFunction likelihood)
{
    return EventSimulator(config).simulate(std::move(likelihood));
}

} // namespace genulens
