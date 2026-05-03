#pragma once

#include "genulens/simulation/run_context.hpp"

namespace genulens {

class Sampler {
public:
    int run_cli(RunContext &context, int argc, char **argv);
};

} // namespace genulens
