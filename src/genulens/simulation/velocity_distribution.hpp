#pragma once

#include <array>
#include "genulens/simulation/run_context.hpp"

namespace genulens {

// Thin interface for velocity sampling. Lifecycle of the underlying Shu DF tables
// (KinematicRuntimeTables::initialize_shu_distribution / release_all) is the
// caller's responsibility.
class VelocityDistribution {
public:
    void bind(RunContext &ctx) { ctx_ = &ctx; }

    // Sample 3D Galactocentric velocity for a star at (D, l, b) with stellar age tau.
    std::array<double, 3> sample(int component, double age,
                                  double D, double l, double b);

private:
    RunContext *ctx_ = nullptr;
};

} // namespace genulens
