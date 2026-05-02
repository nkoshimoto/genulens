#include "genulens/simulation/simulator.hpp"

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

    RandomEngine rng(config_.seed);
    GalacticModel model(config_);
    SimulationResult result;
    const long n = std::max<long>(0, config_.n_simu);
    result.events.reserve(static_cast<std::size_t>(n));

    for (long i = 0; i < n; ++i) {
        const double ds = 4000.0 + 8000.0 * rng.uniform();
        const double dl = ds * rng.uniform();
        const double mass = 0.01 + 1.2 * rng.uniform();
        const double mu = 1.0 + 12.0 * rng.uniform();
        Event event;
        event.tE = std::sqrt(mass) * std::max(1.0, dl * (ds - dl) / ds) / (120.0 * mu);
        event.thetaE = std::sqrt(mass) * 0.3;
        event.piE = event.thetaE > 0.0 ? 0.1 / event.thetaE : 0.0;
        event.lens_distance_pc = dl;
        event.source_distance_pc = ds;
        event.lens_mass_msun = mass;
        event.mu_rel_masyr = mu;
        event.lens_component = static_cast<int>(rng.uniform() * 8.0);
        event.source_component = 8;
        event.weight = likelihood(event) * model.density_at(dl).total();
        result.events.push_back(event);
    }
    return result;
}

SimulationResult simulate(const GenulensConfig &config, LikelihoodFunction likelihood)
{
    return EventSimulator(config).simulate(std::move(likelihood));
}

} // namespace genulens

