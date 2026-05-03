#pragma once

#include <functional>
#include "genulens/model/extinction.hpp"
#include "genulens/simulation/los_density_grid.hpp"
#include "genulens/simulation/mass_function.hpp"
#include "genulens/simulation/observation_config.hpp"
#include "genulens/simulation/run_context.hpp"

namespace genulens {

struct PopulationRuntime;

class EventSampler {
public:
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
    };

    // Run the Monte Carlo loop and print to stdout (CLI mode).
    // Custom likelihood multiplies each event's weight if non-null.
    int run_cli(RunContext &ctx,
                LineOfSightDensityGrid &grid,
                PopulationRuntime &pop,
                MassFunction &mf,
                const Config &cfg,
                const ObservationConfig &obs,
                std::function<double(double /*tE*/, double /*thetaE*/,
                                     double /*piE*/, double /*murel*/,
                                     double /*M_l*/, double /*D_l*/,
                                     double /*D_s*/)> custom_likelihood = nullptr);
};

} // namespace genulens
