#pragma once

#include <functional>
#include <string>
#include <vector>
#include "genulens/model/extinction.hpp"
#include "genulens/model/forward_source.hpp"
#include "genulens/simulation/likelihood.hpp"
#include "genulens/simulation/los_density_grid.hpp"
#include "genulens/simulation/mass_function.hpp"
#include "genulens/simulation/observation_config.hpp"
#include "genulens/simulation/run_context.hpp"

namespace genulens {

struct PopulationRuntime;

class EventSampler {
public:
    using EventSink = std::function<void(const Event &)>;

    struct Config {
        // Precomputed sightline geometry
        double cosPA = 1.0, sinPA = 0.0;
        double cosb  = 1.0, sinb  = 0.0;
        double cosl  = 1.0, sinl  = 0.0;
        double l = 0.0, b = 0.0;

        // Extinction
        const model::ExponentialDustExtinction *extinction = nullptr;
        double AI0 = 0.0, AK0 = 0.0;

        // BH kick
        int    BHhd    = 0,     BHhb  = 0;
        double vkickBH = 100.0, vkickNS = 350.0;
        int    MXDkick = 0;
        double betaBH  = 0.820;

        // Source density (from grid + optional CALCTAU)
        double nallS  = 0.0;
        double Nsall  = 0.0;
        double tauall = 0.0;

        // Reference distance for musRCG correction
        double Dmean = 0.0;

        // Optional source-star sampling from the shared isochrone lookup.
        // Non-empty source_selection_bands are already folded into the
        // line-of-sight source density grid; the event loop samples the
        // concrete source star consistent with that selected source prior.
        const model::ForwardSourceGenerator *forward_source_generator = nullptr;
        RandomEngine *forward_source_rng = nullptr;
        double source_min_initial_mass_msun = 0.09;
        double source_max_initial_mass_msun = 1.0;
        double source_selection_distance_bin_pc = 0.0;
        std::vector<std::string> source_selection_bands;
        std::vector<double> source_selection_min_magnitudes;
        std::vector<double> source_selection_max_magnitudes;
        std::vector<int> source_selection_apparent_magnitudes;
    };

    // Run the Monte Carlo loop. CLI output is controlled by emit_cli_output.
    // Custom likelihood multiplies each event's weight if non-null.
    int run(RunContext &ctx,
            LineOfSightDensityGrid &grid,
            PopulationRuntime &pop,
            MassFunction &mf,
            const Config &cfg,
            const ObservationConfig &obs,
            LikelihoodFunction custom_likelihood = nullptr,
            EventSink event_sink = nullptr,
            bool emit_cli_output = true);

    int run_cli(RunContext &ctx,
                LineOfSightDensityGrid &grid,
                PopulationRuntime &pop,
                MassFunction &mf,
                const Config &cfg,
                const ObservationConfig &obs);
};

} // namespace genulens
