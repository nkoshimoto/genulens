#include "genulens/model/forward_source.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace genulens::model {
namespace {

double interpolate_cumulative(const MassFunctionGrid &grid, double log_mass)
{
    if (log_mass <= grid.log_mass.front()) return grid.cumulative_number_norm.front();
    if (log_mass >= grid.log_mass.back()) return grid.cumulative_number_norm.back();
    const auto upper = std::lower_bound(grid.log_mass.begin(), grid.log_mass.end(), log_mass);
    const std::size_t hi = static_cast<std::size_t>(upper - grid.log_mass.begin());
    const std::size_t lo = hi - 1;
    const double t = (log_mass - grid.log_mass[lo]) / (grid.log_mass[hi] - grid.log_mass[lo]);
    return grid.cumulative_number_norm[lo] +
           (grid.cumulative_number_norm[hi] - grid.cumulative_number_norm[lo]) * t;
}

double invert_cumulative(const MassFunctionGrid &grid, double cumulative)
{
    if (cumulative <= grid.cumulative_number_norm.front()) return grid.mass.front();
    if (cumulative >= grid.cumulative_number_norm.back()) return grid.mass.back();
    const auto upper = std::lower_bound(grid.cumulative_number_norm.begin(),
                                        grid.cumulative_number_norm.end(),
                                        cumulative);
    const std::size_t hi = static_cast<std::size_t>(upper - grid.cumulative_number_norm.begin());
    const std::size_t lo = hi - 1;
    const double t = (cumulative - grid.cumulative_number_norm[lo]) /
                     (grid.cumulative_number_norm[hi] - grid.cumulative_number_norm[lo]);
    const double log_mass = grid.log_mass[lo] + (grid.log_mass[hi] - grid.log_mass[lo]) * t;
    return std::pow(10.0, log_mass);
}

} // namespace

ForwardSourceGenerator::ForwardSourceGenerator(StellarPopulationModel population_model,
                                               IMFParameters imf_parameters)
    : population_model_(std::move(population_model)),
      mass_grid_(BrokenPowerLawIMF(imf_parameters).build_grid(1000,
                                                              imf_parameters.ml,
                                                              imf_parameters.mu))
{
}

ForwardSourceGenerator ForwardSourceGenerator::load_default_roman(IMFParameters imf_parameters)
{
    return ForwardSourceGenerator(StellarPopulationModel::load_default_roman(), imf_parameters);
}

ForwardSourceGenerator ForwardSourceGenerator::load_default_prime(IMFParameters imf_parameters)
{
    return ForwardSourceGenerator(StellarPopulationModel::load_default_prime(), imf_parameters);
}

ForwardSource ForwardSourceGenerator::sample(const ForwardSourceQuery &query, genulens::RandomEngine &rng) const
{
    const double initial_mass = sample_initial_mass(query.min_initial_mass_msun,
                                                   query.max_initial_mass_msun,
                                                   rng);
    StellarPopulationQuery population_query;
    population_query.component = query.component;
    population_query.component_index = query.component_index;
    population_query.initial_mass_msun = initial_mass;
    population_query.log_age = query.log_age;
    population_query.metallicity_mh = query.metallicity_mh;
    population_query.use_default_log_age = query.use_default_log_age;
    population_query.use_default_metallicity = query.use_default_metallicity;

    ForwardSource source;
    source.stellar = population_model_.lookup(population_query);
    source.distance_pc = query.distance_pc;
    source.angular_radius_microarcsec =
        angular_radius_microarcsec(source.stellar.radius_rsun, query.distance_pc);
    return source;
}

ForwardSourceResult ForwardSourceGenerator::sample_many(const ForwardSourceQuery &query,
                                                        std::size_t n_sources,
                                                        genulens::RandomEngine &rng) const
{
    ForwardSourceResult result;
    result.bands = population_model_.isochrones().bands();
    result.sources.reserve(n_sources);
    for (std::size_t i = 0; i < n_sources; ++i) {
        result.sources.push_back(sample(query, rng));
    }
    return result;
}

double ForwardSourceGenerator::sample_initial_mass(double min_mass_msun,
                                                  double max_mass_msun,
                                                  genulens::RandomEngine &rng) const
{
    if (!(min_mass_msun > 0.0) || !(max_mass_msun > min_mass_msun)) {
        throw std::runtime_error("invalid forward-source initial-mass range");
    }
    const double lo = std::max(min_mass_msun, mass_grid_.mass.front());
    const double hi = std::min(max_mass_msun, mass_grid_.mass.back());
    if (!(hi > lo)) {
        throw std::runtime_error("forward-source initial-mass range outside IMF grid");
    }
    const double cmin = interpolate_cumulative(mass_grid_, std::log10(lo));
    const double cmax = interpolate_cumulative(mass_grid_, std::log10(hi));
    const double draw = cmin + (cmax - cmin) * rng.uniform();
    return invert_cumulative(mass_grid_, draw);
}

double angular_radius_microarcsec(double radius_rsun, double distance_pc)
{
    if (!(distance_pc > 0.0)) {
        throw std::runtime_error("angular radius requires positive distance");
    }
    constexpr double kSolarRadiusAu = 0.004650467260962157;
    constexpr double kArcsecToMicroarcsec = 1.0e6;
    return radius_rsun * kSolarRadiusAu / distance_pc * kArcsecToMicroarcsec;
}

std::vector<std::string> ForwardSourceResult::columns() const
{
    std::vector<std::string> out = {
        "iS",
        "D_S",
        "logage_S",
        "MH_S",
        "Zini_S",
        "M_S_ini",
        "M_S",
        "R_S",
        "teff_S",
        "logg_S",
        "theta_S",
    };
    out.reserve(out.size() + bands.size());
    for (const auto &band : bands) {
        out.push_back("M_" + band + "_S");
    }
    return out;
}

std::vector<double> ForwardSourceResult::flattened_rows() const
{
    std::vector<double> rows;
    const auto cols = columns();
    rows.reserve(sources.size() * cols.size());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto &source : sources) {
        const auto &stellar = source.stellar;
        rows.push_back(static_cast<double>(stellar.component_index));
        rows.push_back(source.distance_pc);
        rows.push_back(stellar.log_age);
        rows.push_back(stellar.metallicity_mh);
        rows.push_back(stellar.zini);
        rows.push_back(stellar.initial_mass_msun);
        rows.push_back(stellar.current_mass_msun);
        rows.push_back(stellar.radius_rsun);
        rows.push_back(stellar.teff_k);
        rows.push_back(stellar.logg);
        rows.push_back(source.angular_radius_microarcsec);
        for (const auto &band : bands) {
            const auto found = stellar.absolute_magnitudes.find(band);
            rows.push_back(found == stellar.absolute_magnitudes.end() ? nan : found->second);
        }
    }
    return rows;
}

} // namespace genulens::model
