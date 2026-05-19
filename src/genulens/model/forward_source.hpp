#pragma once

#include "genulens/model/mass_function.hpp"
#include "genulens/model/stellar_population_model.hpp"
#include "genulens/rng.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

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
    std::vector<MagnitudeSelection> magnitude_selections;
};

struct ForwardSource {
    StellarProperties stellar;
    double distance_pc = 0.0;
    double angular_radius_microarcsec = 0.0;
};

struct ForwardSourceResult {
    std::vector<ForwardSource> sources;
    std::vector<std::string> bands;

    std::vector<std::string> columns() const;
    std::vector<double> flattened_rows() const;
    std::size_t row_count() const { return sources.size(); }
    std::size_t column_count() const { return columns().size(); }
};

class ForwardSourceGenerator {
public:
    struct PopulationComponent {
        StellarPopulationModel population_model;
        double weight = 1.0;
        std::string label;
    };

    ForwardSourceGenerator(StellarPopulationModel population_model, IMFParameters imf_parameters);
    ForwardSourceGenerator(std::vector<PopulationComponent> population_components,
                           IMFParameters imf_parameters);

    static ForwardSourceGenerator load_default_roman(IMFParameters imf_parameters = default_model_parameters().imf);
    static ForwardSourceGenerator load_default_prime(IMFParameters imf_parameters = default_model_parameters().imf);
    static ForwardSourceGenerator load_roman(const std::string &primary_table_path,
                                             IMFParameters imf_parameters = default_model_parameters().imf);
    static ForwardSourceGenerator load_prime(const std::string &primary_table_path,
                                             IMFParameters imf_parameters = default_model_parameters().imf);
    static ForwardSourceGenerator load_mixture(const std::string &primary_table_path,
                                               const std::string &secondary_table_path,
                                               double secondary_fraction,
                                               std::vector<double> secondary_fraction_by_component = {},
                                               IMFParameters imf_parameters = default_model_parameters().imf);

    const std::vector<std::string> &bands() const { return population_components_.front().population_model.isochrones().bands(); }
    const std::string &cache_key() const { return cache_key_; }
    ForwardSource sample(const ForwardSourceQuery &query, genulens::RandomEngine &rng) const;
    ForwardSourceResult sample_many(const ForwardSourceQuery &query, std::size_t n_sources,
                                    genulens::RandomEngine &rng) const;
    double selection_probability(const ForwardSourceQuery &query) const;

private:
    struct IntervalCacheState;

    std::vector<PopulationComponent> population_components_;
    std::vector<double> secondary_fraction_by_component_;
    MassFunctionGrid mass_grid_;
    std::string cache_key_;
    std::shared_ptr<IntervalCacheState> interval_cache_;

    std::vector<double> component_weights_for_query(const ForwardSourceQuery &query) const;
    const PopulationComponent &sample_population_component(
        const ForwardSourceQuery &query,
        const std::vector<MagnitudeSelection> &selection,
        genulens::RandomEngine &rng) const;
    ForwardSource sample_from_population(const PopulationComponent &component,
                                         const ForwardSourceQuery &query,
                                         const std::vector<MagnitudeSelection> &selection,
                                         genulens::RandomEngine &rng) const;
    double selection_probability(const PopulationComponent &component,
                                 const ForwardSourceQuery &query) const;
    std::vector<MassInterval> cached_matching_initial_mass_intervals(
        const PopulationComponent &component,
        const ForwardSourceQuery &query) const;
    double sample_initial_mass(double min_mass_msun, double max_mass_msun, genulens::RandomEngine &rng) const;
    double sample_initial_mass_from_interval_bounds(const std::vector<MassInterval> &intervals,
                                                    double min_mass_msun,
                                                    double max_mass_msun,
                                                    genulens::RandomEngine &rng) const;
    double sample_initial_mass_from_intervals(const std::vector<MassInterval> &intervals,
                                              double min_mass_msun,
                                              double max_mass_msun,
                                              genulens::RandomEngine &rng) const;
    double imf_integral(double min_mass_msun, double max_mass_msun) const;
    double imf_integral(const std::vector<MassInterval> &intervals,
                        double min_mass_msun,
                        double max_mass_msun) const;
};

double angular_radius_microarcsec(double radius_rsun, double distance_pc);

} // namespace genulens::model
