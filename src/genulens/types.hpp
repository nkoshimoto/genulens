#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace genulens {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct GalacticCoordinates {
    double l_deg = 1.0;
    double b_deg = -3.9;
};

struct ComponentDensities {
    std::array<double, 11> values{};

    double total() const;
};

struct Event {
    static double missing_value() { return std::numeric_limits<double>::quiet_NaN(); }

    double weight = 1.0;
    double tE = 0.0;
    double thetaE = 0.0;
    double piE = 0.0;
    double lens_distance_pc = 0.0;
    double source_distance_pc = 0.0;
    double lens_mass_msun = 0.0;
    double mu_rel_masyr = 0.0;
    int lens_component = -1;
    int source_component = -1;
    double source_log_age = missing_value();
    double source_metallicity_mh = missing_value();
    double source_zini = missing_value();
    double source_initial_mass_msun = missing_value();
    double source_current_mass_msun = missing_value();
    double source_radius_rsun = missing_value();
    double source_teff_k = missing_value();
    double source_logg = missing_value();
    double source_angular_radius_microarcsec = missing_value();
    std::vector<double> source_absolute_magnitudes;
};

struct SimulationResult {
    std::vector<Event> events;
    bool include_source_properties = false;
    std::vector<std::string> source_property_bands;

    std::vector<std::string> columns() const;
    std::vector<double> flattened_rows() const;
    std::size_t row_count() const { return events.size(); }
    std::size_t column_count() const { return columns().size(); }
};

} // namespace genulens
