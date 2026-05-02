#include "genulens/types.hpp"

#include <numeric>

namespace genulens {

double ComponentDensities::total() const
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

std::vector<std::string> SimulationResult::columns() const
{
    return {
        "weight",
        "tE",
        "thetaE",
        "piE",
        "lens_distance_pc",
        "source_distance_pc",
        "lens_mass_msun",
        "mu_rel_masyr",
        "lens_component",
        "source_component",
    };
}

std::vector<double> SimulationResult::flattened_rows() const
{
    std::vector<double> rows;
    rows.reserve(events.size() * column_count());
    for (const auto &event : events) {
        rows.push_back(event.weight);
        rows.push_back(event.tE);
        rows.push_back(event.thetaE);
        rows.push_back(event.piE);
        rows.push_back(event.lens_distance_pc);
        rows.push_back(event.source_distance_pc);
        rows.push_back(event.lens_mass_msun);
        rows.push_back(event.mu_rel_masyr);
        rows.push_back(static_cast<double>(event.lens_component));
        rows.push_back(static_cast<double>(event.source_component));
    }
    return rows;
}

} // namespace genulens

