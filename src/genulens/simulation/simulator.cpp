#include "genulens/simulation/simulator.hpp"

#include "genulens/simulation/backend.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace genulens {
namespace {

GenulensConfig normalize_source_mode(GenulensConfig config)
{
    if (config.source.mode == "classic") {
        if (config.forward_source.enabled != 0) {
            return config;
        }
        config.forward_source.enabled = 0;
        config.forward_source.selection_bands.clear();
        config.forward_source.selection_min_magnitudes.clear();
        config.forward_source.selection_max_magnitudes.clear();
        config.forward_source.selection_apparent_magnitudes.clear();
        return config;
    }

    if (config.source.mode == "isochrone") {
        if (!(config.source.max_magnitude > config.source.min_magnitude)) {
            throw std::runtime_error("source.max_magnitude must be larger than source.min_magnitude");
        }
        config.source.i_min = config.source.min_magnitude;
        config.source.i_max = config.source.max_magnitude;

        config.forward_source.enabled = 1;
        config.forward_source.photometry = config.source.photometry;
        config.forward_source.isochrone_model = config.source.isochrone_model;
        config.forward_source.isochrone_family = config.source.isochrone_family;
        config.forward_source.isochrone_abundance = config.source.isochrone_abundance;
        config.forward_source.isochrone_alpha_fe = config.source.isochrone_alpha_fe;
        config.forward_source.isochrone_table_path = config.source.isochrone_table_path;
        config.forward_source.secondary_isochrone_family = config.source.secondary_isochrone_family;
        config.forward_source.secondary_isochrone_abundance = config.source.secondary_isochrone_abundance;
        config.forward_source.secondary_isochrone_alpha_fe = config.source.secondary_isochrone_alpha_fe;
        config.forward_source.secondary_isochrone_table_path = config.source.secondary_isochrone_table_path;
        config.forward_source.alpha_enhanced_fraction = config.source.alpha_enhanced_fraction;
        config.forward_source.alpha_enhanced_components = config.source.alpha_enhanced_components;
        config.forward_source.alpha_enhanced_component_fractions =
            config.source.alpha_enhanced_component_fractions;
        config.forward_source.min_initial_mass_msun = config.source.min_initial_mass_msun;
        config.forward_source.max_initial_mass_msun = config.source.max_initial_mass_msun;
        if (config.source.use_magnitude_selection != 0) {
            config.forward_source.selection_bands = {config.source.band};
            config.forward_source.selection_min_magnitudes = {config.source.min_magnitude};
            config.forward_source.selection_max_magnitudes = {config.source.max_magnitude};
            config.forward_source.selection_apparent_magnitudes = {config.source.apparent_magnitude};
        } else {
            config.forward_source.selection_bands.clear();
            config.forward_source.selection_min_magnitudes.clear();
            config.forward_source.selection_max_magnitudes.clear();
            config.forward_source.selection_apparent_magnitudes.clear();
        }
        return config;
    }

    throw std::runtime_error("source.mode must be 'classic' or 'isochrone'");
}

} // namespace

EventSimulator::EventSimulator(GenulensConfig config)
    : config_(normalize_source_mode(std::move(config)))
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
