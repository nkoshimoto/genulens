#pragma once

#include "genulens/rng.hpp"

namespace genulens {

enum class DetectionMode {
    Detection = 0,
    UpperLimit = 1,
    LowerLimit = 2,
};

struct ObservationConstraint {
    double observed = 0.0;
    double error = 0.0;
    double error_inflation = 0.0;
    DetectionMode detection_mode = DetectionMode::Detection;
    bool uniform = false;
};

class ObservationLikelihood {
public:
    explicit ObservationLikelihood(ObservationConstraint constraint);

    double accept_weight(double model_value, RandomEngine &rng) const;

private:
    ObservationConstraint constraint_;
};

} // namespace genulens

