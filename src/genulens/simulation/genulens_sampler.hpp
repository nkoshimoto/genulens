#pragma once

#include "genulens/simulation/run_context.hpp"

namespace genulens {

class GenulensSampler {
public:
    int run_cli(GenulensRunContext &context, int argc, char **argv);
};

} // namespace genulens
