#pragma once

#include "genulens/simulation/scientific_state.hpp"

namespace genulens {

class ScientificEngine {
public:
    int run(int argc, char **argv);

private:
    ScientificState state_;
};

} // namespace genulens

