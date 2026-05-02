#pragma once

#include "genulens/model/density.hpp"
#include "genulens/model/kinematics.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/options.hpp"

namespace genulens {

class GalacticModel {
public:
    explicit GalacticModel(GenulensConfig config = {});

    const GenulensConfig &config() const { return config_; }
    GalacticCoordinates coordinates() const { return {config_.l, config_.b}; }
    ComponentDensities density_at(double distance_pc) const;
    Vec3 cartesian_at(double distance_pc) const;
    double imf(double mass_msun) const;

private:
    GenulensConfig config_;
};

} // namespace genulens

