#pragma once

#include "genulens/model/isochrone_grid.hpp"

#include <string>

namespace genulens::model {

struct StellarPopulationQuery {
    std::string component;
    int component_index = -1;
    double initial_mass_msun = 0.0;
    double log_age = 0.0;
    double metallicity_mh = 0.0;
    bool use_default_log_age = true;
    bool use_default_metallicity = true;
};

class StellarPopulationModel {
public:
    static StellarPopulationModel load_default_roman();
    static StellarPopulationModel load_default_prime();

    explicit StellarPopulationModel(IsochroneGrid grid);

    const IsochroneGrid &isochrones() const { return isochrones_; }
    StellarProperties lookup(const StellarPopulationQuery &query) const;

    static std::string component_name(int component_index);
    static int component_index(const std::string &component);
    static double default_log_age(const std::string &component);
    static double default_metallicity_mh(const std::string &component);

private:
    IsochroneGrid isochrones_;
};

} // namespace genulens::model
