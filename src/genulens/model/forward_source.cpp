#include "genulens/model/forward_source.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace genulens::model {
struct ForwardSourceGenerator::IntervalCacheState {
    std::mutex mutex;
    std::unordered_map<std::string, std::vector<MassInterval>> intervals;
};

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

bool matches_selection(const StellarProperties &stellar,
                       const std::vector<MagnitudeSelection> &selection)
{
    for (const auto &cut : selection) {
        const auto found = stellar.absolute_magnitudes.find(cut.band);
        if (found == stellar.absolute_magnitudes.end()) return false;
        const double mag = found->second + cut.magnitude_offset;
        if (!(mag >= cut.min_magnitude && mag <= cut.max_magnitude)) return false;
    }
    return true;
}

std::string interval_cache_key(const ForwardSourceGenerator::PopulationComponent &component,
                               const ForwardSourceQuery &query)
{
    std::ostringstream out;
    out << std::setprecision(17)
        << component.label
        << "|component=" << query.component
        << "|component_index=" << query.component_index
        << "|log_age=" << query.log_age
        << "|mh=" << query.metallicity_mh
        << "|default_age=" << query.use_default_log_age
        << "|default_mh=" << query.use_default_metallicity
        << "|selection=";
    for (const auto &selection : query.magnitude_selections) {
        out << selection.band << ','
            << selection.min_magnitude << ','
            << selection.max_magnitude << ','
            << selection.magnitude_offset << ';';
    }
    return out.str();
}

std::string imf_cache_key(const IMFParameters &imf)
{
    std::ostringstream out;
    out << std::setprecision(17)
        << "imf:"
        << imf.m0 << ','
        << imf.m1 << ','
        << imf.m2 << ','
        << imf.m3 << ','
        << imf.ml << ','
        << imf.mu << ','
        << imf.alpha0 << ','
        << imf.alpha1 << ','
        << imf.alpha2 << ','
        << imf.alpha3 << ','
        << imf.alpha4;
    return out.str();
}

} // namespace

ForwardSourceGenerator::ForwardSourceGenerator(StellarPopulationModel population_model,
                                               IMFParameters imf_parameters)
    : ForwardSourceGenerator(
          std::vector<PopulationComponent>{{std::move(population_model), 1.0, "primary"}},
          imf_parameters)
{
}

ForwardSourceGenerator::ForwardSourceGenerator(std::vector<PopulationComponent> population_components,
                                               IMFParameters imf_parameters)
    : population_components_(std::move(population_components)),
      mass_grid_(BrokenPowerLawIMF(imf_parameters).build_grid(1000,
                                                              imf_parameters.ml,
                                                              imf_parameters.mu)),
      cache_key_("anonymous|" + imf_cache_key(imf_parameters)),
      interval_cache_(std::make_shared<IntervalCacheState>())
{
    if (population_components_.empty()) {
        throw std::runtime_error("forward-source generator requires at least one isochrone population");
    }
    double total_weight = 0.0;
    const auto &reference_bands = population_components_.front().population_model.isochrones().bands();
    for (const auto &component : population_components_) {
        if (!(component.weight >= 0.0)) {
            throw std::runtime_error("forward-source isochrone population weight must be non-negative");
        }
        total_weight += component.weight;
        if (component.population_model.isochrones().bands() != reference_bands) {
            throw std::runtime_error("mixed forward-source isochrone populations must use the same bands");
        }
    }
    if (!(total_weight > 0.0)) {
        throw std::runtime_error("forward-source isochrone population weights sum to zero");
    }
    for (auto &component : population_components_) {
        component.weight /= total_weight;
    }
}

ForwardSourceGenerator ForwardSourceGenerator::load_default_roman(IMFParameters imf_parameters)
{
    auto generator = ForwardSourceGenerator(StellarPopulationModel::load_default_roman(), imf_parameters);
    generator.cache_key_ = "default_roman|" + imf_cache_key(imf_parameters);
    return generator;
}

ForwardSourceGenerator ForwardSourceGenerator::load_default_prime(IMFParameters imf_parameters)
{
    auto generator = ForwardSourceGenerator(StellarPopulationModel::load_default_prime(), imf_parameters);
    generator.cache_key_ = "default_prime|" + imf_cache_key(imf_parameters);
    return generator;
}

ForwardSourceGenerator ForwardSourceGenerator::load_roman(const std::string &primary_table_path,
                                                          IMFParameters imf_parameters)
{
    auto generator = ForwardSourceGenerator(
        StellarPopulationModel(IsochroneGrid::load(primary_table_path)), imf_parameters);
    generator.cache_key_ = "roman|" + primary_table_path + "|" + imf_cache_key(imf_parameters);
    return generator;
}

ForwardSourceGenerator ForwardSourceGenerator::load_prime(const std::string &primary_table_path,
                                                          IMFParameters imf_parameters)
{
    auto generator = ForwardSourceGenerator(
        StellarPopulationModel(IsochroneGrid::load(primary_table_path)), imf_parameters);
    generator.cache_key_ = "prime|" + primary_table_path + "|" + imf_cache_key(imf_parameters);
    return generator;
}

ForwardSourceGenerator ForwardSourceGenerator::load_mixture(const std::string &primary_table_path,
                                                            const std::string &secondary_table_path,
                                                            double secondary_fraction,
                                                            std::vector<double> secondary_fraction_by_component,
                                                            IMFParameters imf_parameters)
{
    if (!(secondary_fraction >= 0.0 && secondary_fraction <= 1.0)) {
        throw std::runtime_error("secondary isochrone fraction must be in [0, 1]");
    }
    std::vector<PopulationComponent> components;
    components.push_back({
        StellarPopulationModel(IsochroneGrid::load(primary_table_path)),
        1.0 - secondary_fraction,
        "primary",
    });
    components.push_back({
        StellarPopulationModel(IsochroneGrid::load(secondary_table_path)),
        secondary_fraction,
        "secondary",
    });
    auto generator = ForwardSourceGenerator(std::move(components), imf_parameters);
    if (!secondary_fraction_by_component.empty()) {
        if (secondary_fraction_by_component.size() < 11) {
            throw std::runtime_error("component alpha-enhanced fraction vector must cover components 0..10");
        }
        for (const double fraction : secondary_fraction_by_component) {
            if (!(fraction >= 0.0 && fraction <= 1.0)) {
                throw std::runtime_error("component alpha-enhanced fractions must be in [0, 1]");
            }
        }
        generator.secondary_fraction_by_component_ = std::move(secondary_fraction_by_component);
    }
    std::ostringstream key;
    key << std::setprecision(17)
        << "mixture|primary=" << primary_table_path
        << "|secondary=" << secondary_table_path
        << "|secondary_fraction=" << secondary_fraction
        << '|' << imf_cache_key(imf_parameters);
    if (!generator.secondary_fraction_by_component_.empty()) {
        key << "|component_fractions=";
        for (const double fraction : generator.secondary_fraction_by_component_) {
            key << fraction << ',';
        }
    }
    generator.cache_key_ = key.str();
    return generator;
}

ForwardSource ForwardSourceGenerator::sample(const ForwardSourceQuery &query, genulens::RandomEngine &rng) const
{
    if (!query.magnitude_selections.empty()) {
        const auto &component = sample_population_component(query, query.magnitude_selections, rng);
        return sample_from_population(component, query, query.magnitude_selections, rng);
    }
    const std::vector<MagnitudeSelection> no_selection;
    const auto &component = sample_population_component(query, no_selection, rng);
    return sample_from_population(component, query, no_selection, rng);
}

ForwardSourceResult ForwardSourceGenerator::sample_many(const ForwardSourceQuery &query,
                                                        std::size_t n_sources,
                                                        genulens::RandomEngine &rng) const
{
    ForwardSourceResult result;
    result.bands = bands();
    result.sources.reserve(n_sources);
    for (std::size_t i = 0; i < n_sources; ++i) {
        result.sources.push_back(sample(query, rng));
    }
    return result;
}

double ForwardSourceGenerator::selection_probability(const ForwardSourceQuery &query) const
{
    if (query.magnitude_selections.empty()) return 1.0;

    double probability = 0.0;
    const auto weights = component_weights_for_query(query);
    for (std::size_t i = 0; i < population_components_.size(); ++i) {
        if (weights[i] <= 0.0) continue;
        probability += weights[i] * selection_probability(population_components_[i], query);
    }
    return std::max(0.0, std::min(1.0, probability));
}

std::vector<double> ForwardSourceGenerator::component_weights_for_query(const ForwardSourceQuery &query) const
{
    std::vector<double> weights;
    weights.reserve(population_components_.size());
    for (const auto &component : population_components_) weights.push_back(component.weight);

    if (population_components_.size() == 2 && !secondary_fraction_by_component_.empty()) {
        int component_index = query.component_index;
        if (component_index == 10) component_index = 7;
        if (component_index >= 0 &&
            static_cast<std::size_t>(component_index) < secondary_fraction_by_component_.size()) {
            const double secondary_fraction =
                secondary_fraction_by_component_[static_cast<std::size_t>(component_index)];
            weights[0] = 1.0 - secondary_fraction;
            weights[1] = secondary_fraction;
        }
    }
    return weights;
}

const ForwardSourceGenerator::PopulationComponent &ForwardSourceGenerator::sample_population_component(
    const ForwardSourceQuery &query,
    const std::vector<MagnitudeSelection> &selection,
    genulens::RandomEngine &rng) const
{
    if (population_components_.size() == 1) return population_components_.front();

    double total = 0.0;
    std::vector<double> cumulative;
    cumulative.reserve(population_components_.size());
    const auto base_weights = component_weights_for_query(query);
    for (std::size_t i = 0; i < population_components_.size(); ++i) {
        const auto &component = population_components_[i];
        double weight = base_weights[i];
        if (weight <= 0.0) {
            cumulative.push_back(total);
            continue;
        }
        if (!selection.empty()) {
            ForwardSourceQuery selected_query = query;
            selected_query.magnitude_selections = selection;
            weight *= selection_probability(component, selected_query);
        }
        total += weight;
        cumulative.push_back(total);
    }
    if (!(total > 0.0)) {
        throw std::runtime_error("forward-source isochrone mixture has zero selected probability");
    }
    const double draw = rng.uniform() * total;
    for (std::size_t i = 0; i < cumulative.size(); ++i) {
        if (draw <= cumulative[i]) return population_components_[i];
    }
    return population_components_.back();
}

ForwardSource ForwardSourceGenerator::sample_from_population(
    const PopulationComponent &component,
    const ForwardSourceQuery &query,
    const std::vector<MagnitudeSelection> &selection,
    genulens::RandomEngine &rng) const
{
    StellarPopulationQuery population_query;
    population_query.component = query.component;
    population_query.component_index = query.component_index;
    population_query.log_age = query.log_age;
    population_query.metallicity_mh = query.metallicity_mh;
    population_query.use_default_log_age = query.use_default_log_age;
    population_query.use_default_metallicity = query.use_default_metallicity;

    std::vector<MassInterval> intervals;
    if (!selection.empty()) {
        intervals = cached_matching_initial_mass_intervals(component, query);
    }
    const int max_attempts = selection.empty() ? 1 : 64;
    for (int attempt = 0; attempt <= max_attempts; ++attempt) {
        double initial_mass = 0.0;
        if (selection.empty()) {
            initial_mass = sample_initial_mass(query.min_initial_mass_msun,
                                               query.max_initial_mass_msun,
                                               rng);
        } else if (attempt < max_attempts) {
            initial_mass = sample_initial_mass_from_interval_bounds(intervals,
                                                                    query.min_initial_mass_msun,
                                                                    query.max_initial_mass_msun,
                                                                    rng);
        } else {
            // Conservative bounds can include gaps around non-monotonic
            // isochrone phases. Fall back to the exact allowed intervals only
            // after bounded proposals fail, without rejecting the event.
            initial_mass = sample_initial_mass_from_intervals(intervals,
                                                              query.min_initial_mass_msun,
                                                              query.max_initial_mass_msun,
                                                              rng);
        }
        population_query.initial_mass_msun = initial_mass;

        ForwardSource source;
        source.stellar = component.population_model.lookup(population_query);
        if (!selection.empty() && !matches_selection(source.stellar, selection)) {
            continue;
        }
        source.distance_pc = query.distance_pc;
        source.angular_radius_microarcsec =
            angular_radius_microarcsec(source.stellar.radius_rsun, query.distance_pc);
        return source;
    }
    throw std::runtime_error("forward-source magnitude selection failed after bounded mass proposals");
}

double ForwardSourceGenerator::selection_probability(const PopulationComponent &component,
                                                     const ForwardSourceQuery &query) const
{
    if (query.magnitude_selections.empty()) return 1.0;

    const auto intervals = cached_matching_initial_mass_intervals(component, query);
    const double total = imf_integral(query.min_initial_mass_msun,
                                      query.max_initial_mass_msun);
    if (!(total > 0.0)) return 0.0;
    const double selected = imf_integral(intervals,
                                         query.min_initial_mass_msun,
                                         query.max_initial_mass_msun);
    return std::max(0.0, std::min(1.0, selected / total));
}

std::vector<MassInterval> ForwardSourceGenerator::cached_matching_initial_mass_intervals(
    const PopulationComponent &component,
    const ForwardSourceQuery &query) const
{
    if (query.magnitude_selections.empty()) return {};
    const auto key = interval_cache_key(component, query);
    {
        const std::lock_guard<std::mutex> lock(interval_cache_->mutex);
        const auto found = interval_cache_->intervals.find(key);
        if (found != interval_cache_->intervals.end()) return found->second;
    }

    StellarPopulationQuery population_query;
    population_query.component = query.component;
    population_query.component_index = query.component_index;
    population_query.log_age = query.log_age;
    population_query.metallicity_mh = query.metallicity_mh;
    population_query.use_default_log_age = query.use_default_log_age;
    population_query.use_default_metallicity = query.use_default_metallicity;

    auto intervals = component.population_model.matching_initial_mass_intervals(
        population_query, query.magnitude_selections);
    {
        const std::lock_guard<std::mutex> lock(interval_cache_->mutex);
        if (interval_cache_->intervals.size() > 200000) interval_cache_->intervals.clear();
        interval_cache_->intervals.emplace(key, intervals);
    }
    return intervals;
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

double ForwardSourceGenerator::imf_integral(double min_mass_msun,
                                            double max_mass_msun) const
{
    if (!(min_mass_msun > 0.0) || !(max_mass_msun > min_mass_msun)) {
        return 0.0;
    }
    const double lo = std::max(min_mass_msun, mass_grid_.mass.front());
    const double hi = std::min(max_mass_msun, mass_grid_.mass.back());
    if (!(hi > lo)) return 0.0;
    const double cmin = interpolate_cumulative(mass_grid_, std::log10(lo));
    const double cmax = interpolate_cumulative(mass_grid_, std::log10(hi));
    return std::max(0.0, cmax - cmin);
}

double ForwardSourceGenerator::imf_integral(const std::vector<MassInterval> &intervals,
                                            double min_mass_msun,
                                            double max_mass_msun) const
{
    double total = 0.0;
    const double lo = std::max(min_mass_msun, mass_grid_.mass.front());
    const double hi = std::min(max_mass_msun, mass_grid_.mass.back());
    if (!(hi > lo)) return 0.0;
    for (const auto &interval : intervals) {
        const double ilo = std::max(lo, interval.min_mass_msun);
        const double ihi = std::min(hi, interval.max_mass_msun);
        if (!(ihi > ilo)) continue;
        total += imf_integral(ilo, ihi);
    }
    return total;
}

double ForwardSourceGenerator::sample_initial_mass_from_interval_bounds(
    const std::vector<MassInterval> &intervals,
    double min_mass_msun,
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

    double bound_lo = std::numeric_limits<double>::infinity();
    double bound_hi = -std::numeric_limits<double>::infinity();
    for (const auto &interval : intervals) {
        const double ilo = std::max(lo, interval.min_mass_msun);
        const double ihi = std::min(hi, interval.max_mass_msun);
        if (!(ihi > ilo)) continue;
        bound_lo = std::min(bound_lo, ilo);
        bound_hi = std::max(bound_hi, ihi);
    }
    if (!(bound_hi > bound_lo)) {
        throw std::runtime_error("forward-source magnitude selection has zero IMF probability");
    }
    return sample_initial_mass(bound_lo, bound_hi, rng);
}

double ForwardSourceGenerator::sample_initial_mass_from_intervals(
    const std::vector<MassInterval> &intervals,
    double min_mass_msun,
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

    struct WeightedInterval {
        double cmin = 0.0;
        double cmax = 0.0;
        double cumulative_weight = 0.0;
    };
    std::vector<WeightedInterval> weighted;
    double allowed = 0.0;
    for (const auto &interval : intervals) {
        const double ilo = std::max(lo, interval.min_mass_msun);
        const double ihi = std::min(hi, interval.max_mass_msun);
        if (!(ihi > ilo)) continue;
        const double icmin = interpolate_cumulative(mass_grid_, std::log10(ilo));
        const double icmax = interpolate_cumulative(mass_grid_, std::log10(ihi));
        if (!(icmax > icmin)) continue;
        allowed += icmax - icmin;
        weighted.push_back({icmin, icmax, allowed});
    }

    if (!(allowed > 0.0)) {
        throw std::runtime_error("forward-source magnitude selection has zero IMF probability");
    }

    const double draw = rng.uniform() * allowed;
    const auto chosen = std::lower_bound(weighted.begin(), weighted.end(), draw,
                                         [](const WeightedInterval &interval, double value) {
                                             return interval.cumulative_weight < value;
                                         });
    const auto &interval = (chosen == weighted.end()) ? weighted.back() : *chosen;
    const double previous = interval.cumulative_weight - (interval.cmax - interval.cmin);
    const double local = interval.cmin + (draw - previous);
    return invert_cumulative(mass_grid_, local);
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
