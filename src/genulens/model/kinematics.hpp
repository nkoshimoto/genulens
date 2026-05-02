#pragma once

#include "genulens/types.hpp"

namespace genulens::model {

Vec3 galactic_to_cartesian(double distance_pc, GalacticCoordinates coordinates, double sun_radius_pc = 8160.0);
double position_angle_deg(double l_deg, double b_deg);

} // namespace genulens::model

