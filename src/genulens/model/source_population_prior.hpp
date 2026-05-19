#pragma once

#include <vector>

namespace genulens::model {

struct AgeMetallicityPoint {
    double log_age = 0.0;
    double metallicity_mh = 0.0;
    double weight = 1.0;
};

class SourcePopulationPrior {
public:
    static std::vector<AgeMetallicityPoint> points_for_component(int component_index);
};

} // namespace genulens::model
