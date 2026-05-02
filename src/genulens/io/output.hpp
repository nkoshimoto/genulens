#pragma once

#include "genulens/types.hpp"

#include <iosfwd>

namespace genulens {

void write_events_tsv(std::ostream &out, const SimulationResult &result);

} // namespace genulens

