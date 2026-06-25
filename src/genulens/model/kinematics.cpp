#include "genulens/model/kinematics.hpp"

#include "genulens/constants.hpp"

#include <cmath>

namespace genulens::model {

Vec3 galactic_to_cartesian(double distance_pc, GalacticCoordinates coordinates, double sun_radius_pc)
{
    const double l = coordinates.l_deg * kPi / 180.0;
    const double b = coordinates.b_deg * kPi / 180.0;
    return {
        sun_radius_pc - distance_pc * std::cos(b) * std::cos(l),
        distance_pc * std::cos(b) * std::sin(l),
        distance_pc * std::sin(b) + kSunZPc,
    };
}

double position_angle_deg(double l_deg, double b_deg)
{
    const double l = l_deg * kPi / 180.0;
    const double b = b_deg * kPi / 180.0;
    return std::atan2(std::sin(l), std::tan(b)) * 180.0 / kPi;
}

} // namespace genulens::model

