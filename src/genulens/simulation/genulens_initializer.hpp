#pragma once

#include "genulens/simulation/run_context.hpp"

namespace genulens {

class GenulensInitializer {
public:
    GenulensRunContext create_context() const;
};

} // namespace genulens
