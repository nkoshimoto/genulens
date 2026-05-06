#pragma once

#include "genulens/model/mass_function.hpp"
#include "genulens/model/stellar_population_model.hpp"
#include "genulens/rng.hpp"

#include <string>

namespace genulens::model {

struct ForwardSourceQuery {
    std::string component;
    int component_index = -1;
    double distance_pc = 8000.0;
    double min_initial_mass_msun = 0.09;
    double max_initial_mass_msun = 1.0;
    double log_age = 0.0;
    double metallicity_mh = 0.0;
    bool use_default_log_age = true;
    bool use_default_metallicity = true;
};

struct ForwardSource {
    StellarProperties stellar;
    double distance_pc = 0.0;
    double angular_radius_microarcsec = 0.0;
};

class ForwardSourceGenerator {
public:
    ForwardSourceGenerator(StellarPopulationModel population_model, IMFParameters imf_parameters);

    static ForwardSourceGenerator load_default_roman(IMFParameters imf_parameters = default_model_parameters().imf);
    static ForwardSourceGenerator load_default_prime(IMFParameters imf_parameters = default_model_parameters().imf);

    ForwardSource sample(const ForwardSourceQuery &query, genulens::RandomEngine &rng) const;

private:
    StellarPopulationModel population_model_;
    MassFunctionGrid mass_grid_;

    double sample_initial_mass(double min_mass_msun, double max_mass_msun, genulens::RandomEngine &rng) const;
};

double angular_radius_microarcsec(double radius_rsun, double distance_pc);

} // namespace genulens::model
