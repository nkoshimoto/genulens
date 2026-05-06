#include "genulens/types.hpp"

#include <numeric>

namespace genulens {

double ComponentDensities::total() const
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

std::vector<std::string> SimulationResult::columns() const
{
    std::vector<std::string> out = {
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
    if (include_source_properties) {
        out.push_back("source_log_age");
        out.push_back("source_metallicity_mh");
        out.push_back("source_zini");
        out.push_back("source_initial_mass_msun");
        out.push_back("source_current_mass_msun");
        out.push_back("source_radius_rsun");
        out.push_back("source_teff_k");
        out.push_back("source_logg");
        out.push_back("source_angular_radius_microarcsec");
        for (const auto &band : source_property_bands) {
            out.push_back("source_abs_" + band);
        }
    }
    return out;
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
        if (include_source_properties) {
            rows.push_back(event.source_log_age);
            rows.push_back(event.source_metallicity_mh);
            rows.push_back(event.source_zini);
            rows.push_back(event.source_initial_mass_msun);
            rows.push_back(event.source_current_mass_msun);
            rows.push_back(event.source_radius_rsun);
            rows.push_back(event.source_teff_k);
            rows.push_back(event.source_logg);
            rows.push_back(event.source_angular_radius_microarcsec);
            for (std::size_t i = 0; i < source_property_bands.size(); ++i) {
                rows.push_back(i < event.source_absolute_magnitudes.size()
                                   ? event.source_absolute_magnitudes[i]
                                   : Event::missing_value());
            }
        }
    }
    return rows;
}

} // namespace genulens
