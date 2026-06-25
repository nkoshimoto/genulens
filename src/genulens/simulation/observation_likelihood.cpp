#include "genulens/simulation/observation_likelihood.hpp"

#include <cmath>

namespace genulens {

ObservationLikelihood::ObservationLikelihood(ObservationConstraint constraint)
    : constraint_(constraint)
{
}

double ObservationLikelihood::accept_weight(double model_value, RandomEngine &rng) const
{
    double gamma = 1.0;
    const bool accepted_by_limit =
        constraint_.detection_mode == DetectionMode::Detection ||
        (constraint_.detection_mode == DetectionMode::UpperLimit && model_value > constraint_.observed) ||
        (constraint_.detection_mode == DetectionMode::LowerLimit && model_value < constraint_.observed);

    if (constraint_.error > 0.0 && accepted_by_limit) {
        if (constraint_.uniform) {
            gamma = (model_value > constraint_.observed - constraint_.error &&
                     model_value < constraint_.observed + constraint_.error)
                        ? 1.0
                        : 0.0;
        } else {
            double chi = (model_value - constraint_.observed) / constraint_.error;
            double reference_like = 1.0;
            if (constraint_.error_inflation > 0.0) {
                reference_like = std::exp(-0.5 * chi * chi);
                chi /= constraint_.error_inflation;
            }
            const double proposal_like = std::exp(-0.5 * chi * chi);
            gamma = (proposal_like > rng.uniform()) ? 1.0 : 0.0;
            if (constraint_.error_inflation > 0.0) {
                gamma *= reference_like / proposal_like;
            }
        }
    } else if (!accepted_by_limit) {
        gamma = 0.0;
    }
    return gamma;
}

} // namespace genulens

