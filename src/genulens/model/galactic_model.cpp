#include "genulens/model/galactic_model.hpp"

namespace genulens {

GalacticModel::GalacticModel(GenulensConfig config)
    : config_(std::move(config))
{
}

ComponentDensities GalacticModel::density_at(double distance_pc) const
{
    return model::approximate_density(distance_pc, coordinates());
}

Vec3 GalacticModel::cartesian_at(double distance_pc) const
{
    return model::galactic_to_cartesian(distance_pc, coordinates());
}

double GalacticModel::imf(double mass_msun) const
{
    return model::broken_power_law_imf(mass_msun, config_.model.imf);
}

} // namespace genulens
