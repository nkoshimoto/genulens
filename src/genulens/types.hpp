#pragma once

#include <array>
#include <cstddef>
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
};

struct SimulationResult {
    std::vector<Event> events;

    std::vector<std::string> columns() const;
    std::vector<double> flattened_rows() const;
    std::size_t row_count() const { return events.size(); }
    std::size_t column_count() const { return columns().size(); }
};

} // namespace genulens

