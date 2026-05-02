#include "genulens/io/output.hpp"

#include <ostream>

namespace genulens {

void write_events_tsv(std::ostream &out, const SimulationResult &result)
{
    const auto names = result.columns();
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (i) out << '\t';
        out << names[i];
    }
    out << '\n';
    for (const auto &event : result.events) {
        out << event.weight << '\t'
            << event.tE << '\t'
            << event.thetaE << '\t'
            << event.piE << '\t'
            << event.lens_distance_pc << '\t'
            << event.source_distance_pc << '\t'
            << event.lens_mass_msun << '\t'
            << event.mu_rel_masyr << '\t'
            << event.lens_component << '\t'
            << event.source_component << '\n';
    }
}

} // namespace genulens

