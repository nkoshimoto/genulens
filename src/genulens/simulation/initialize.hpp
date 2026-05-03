#pragma once

#include "genulens/simulation/run_context.hpp"

namespace genulens {

class Initializer {
public:
    RunContext create_context() const;
    void initialize_rng(RunContext &context, int argc, char **argv) const;
    void read_model_options(RunContext &context, int argc, char **argv) const;
    void finalize_spatial_model(RunContext &context, int argc, char **argv) const;
    void read_sampling_options(RunContext &context, int argc, char **argv,
                               double cos_pa, double sin_pa) const;
};

} // namespace genulens
