#pragma once

#include <array>

namespace genulens::model {

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

