#include "genulens/model/extinction.hpp"

#include <cmath>

namespace genulens::model {

ExponentialDustExtinction::ExponentialDustExtinction(double scale_pc, double reference_distance_pc,
                                                     double ai_reference, double ak_reference,
                                                     double evi_reference)
    : scale_pc_(scale_pc),
      reference_distance_pc_(reference_distance_pc)
{
    const double denominator = (scale_pc_ > 0.0 && reference_distance_pc_ > 0.0)
                                   ? 1.0 - std::exp(-reference_distance_pc_ / scale_pc_)
                                   : 0.0;
    ai0_ = denominator != 0.0 ? ai_reference / denominator : 0.0;
    ak0_ = denominator != 0.0 ? ak_reference / denominator : 0.0;
    evi0_ = denominator != 0.0 ? evi_reference / denominator : 0.0;
}

BandExtinction ExponentialDustExtinction::at_distance(double distance_pc) const
{
    const double factor = 1.0 - std::exp(-distance_pc / scale_pc_);
    return {ai0_ * factor, ak0_ * factor, evi0_ * factor};
}

double ExponentialDustExtinction::distance_modulus_term(double distance_pc) const
{
    return 5.0 * std::log10(0.1 * distance_pc);
}

double source_distance_weight(double distance_pc, double gamma_ds)
{
    return distance_pc > 0.0 ? std::pow(distance_pc / 8000.0, gamma_ds) : 0.0;
}

} // namespace genulens::model
