#include "genulens/model/extinction.hpp"

#include <cmath>

namespace genulens::model {

double source_distance_weight(double distance_pc, double gamma_ds)
{
    return distance_pc > 0.0 ? std::pow(distance_pc / 8000.0, gamma_ds) : 0.0;
}

} // namespace genulens::model

