#pragma once

#include "genulens/simulation/event_sampler.hpp"
#include "genulens/simulation/los_density_grid.hpp"
#include "genulens/simulation/mass_function.hpp"
#include "genulens/simulation/observation_config.hpp"
#include "genulens/simulation/internal/runtime.hpp"

#include <memory>

namespace genulens {

struct PreparedSimulation {
    PopulationRuntime population;
    KinematicRuntimeTables kinematic_tables;
    NsdMomentRuntime nsd_moments;
    LineOfSightDensityGrid grid;
    MassFunction mass_function;
    std::unique_ptr<model::ForwardSourceGenerator> forward_source_generator;
    std::unique_ptr<RandomEngine> forward_source_rng;
    EventSampler::Config event_config;
    ObservationConfig observation;

    bool luminosity_functions_active = false;
    bool population_active = false;
    bool kinematic_tables_active = false;
    bool nsd_moments_active = false;
    bool density_sightline_allocated = false;

    void cleanup(RunContext &context);
};

int run_prepared_events(RunContext &context,
                        PreparedSimulation &prepared,
                        LikelihoodFunction custom_likelihood,
                        EventSampler::EventSink event_sink,
                        bool emit_cli_output);

} // namespace genulens
