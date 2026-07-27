#pragma once

#include <array>

namespace genulens::model {

// Return Galactic longitude in the conventional [-180, 180) interval.
// The Galactic model accepts either 0--360 or signed longitudes, but
// empirical relations written around the Galactic centre require the
// signed representation.
double wrapped_longitude_deg(double l_deg);

struct PositionAngle {
    double degrees = 0.0;
    double cos_pa = 1.0;
    double sin_pa = 0.0;
};

class CoordinateTransformer {
public:
    explicit CoordinateTransformer(std::array<double, 3> sgr_a_offset = {});

    std::array<double, 3> distance_l_b_to_xyz(double distance_pc, double l_deg, double b_deg,
                                              double sun_radius_pc) const;
    static PositionAngle position_angle(double l_deg, double b_deg);

private:
    std::array<double, 3> sgr_a_offset_{};
};

} // namespace genulens::model
