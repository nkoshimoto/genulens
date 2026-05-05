#include "genulens/simulation/prepared_simulation.hpp"

#include <cstdlib>
#include <utility>

namespace genulens {

void PreparedSimulation::cleanup(RunContext &context)
{
    if (nsd_moments_active) {
        nsd_moments.release_if_enabled(context);
        nsd_moments_active = false;
    }
    if (luminosity_functions_active) {
        population.release_luminosity_functions(context);
        luminosity_functions_active = false;
    }
    if (population_active) {
        population.release_all(context);
        population_active = false;
    }
    if (density_sightline_allocated) {
        free(context.density.lDs);
        free(context.density.bDs);
        context.density.lDs = nullptr;
        context.density.bDs = nullptr;
        density_sightline_allocated = false;
    }
    if (kinematic_tables_active) {
        kinematic_tables.release_all(context);
        kinematic_tables_active = false;
    }
}

int run_prepared_events(RunContext &context,
                        PreparedSimulation &prepared,
                        LikelihoodFunction custom_likelihood,
                        EventSampler::EventSink event_sink,
                        bool emit_cli_output)
{
    EventSampler event_sampler;
    return event_sampler.run(context,
                             prepared.grid,
                             prepared.population,
                             prepared.mass_function,
                             prepared.event_config,
                             prepared.observation,
                             std::move(custom_likelihood),
                             std::move(event_sink),
                             emit_cli_output);
}

} // namespace genulens
