#include "genulens/model/stellar_population_model.hpp"

#include <array>
#include <stdexcept>
#include <utility>

namespace genulens::model {
namespace {

struct ComponentDefaults {
    const char *name;
    int index;
    double log_age;
    double metallicity_mh;
};

constexpr std::array<ComponentDefaults, 10> kDefaults = {{
    {"thin1", 0, 8.0, 0.0},
    {"thin2", 1, 8.77815125, 0.0},
    {"thin3", 2, 9.176091259, 0.0},
    {"thin4", 3, 9.397940009, 0.0},
    {"thin5", 4, 9.607455023, 0.0},
    {"thin6", 5, 9.781755375, 0.0},
    {"thin7", 6, 9.937016107, 0.0},
    {"thick", 7, 10.079181246047625, -0.8},
    {"bar", 8, 9.95424, 0.0},
    {"NSD", 9, 9.8451, 0.0},
}};

const ComponentDefaults &find_component(const std::string &component)
{
    for (const auto &entry : kDefaults) {
        if (component == entry.name) return entry;
    }
    throw std::runtime_error("unknown stellar component: " + component);
}

const ComponentDefaults &find_component(int component_index)
{
    for (const auto &entry : kDefaults) {
        if (component_index == entry.index) return entry;
    }
    throw std::runtime_error("unknown stellar component index: " + std::to_string(component_index));
}

} // namespace

StellarPopulationModel StellarPopulationModel::load_default_roman()
{
    return StellarPopulationModel(IsochroneGrid::load_default_roman());
}

StellarPopulationModel StellarPopulationModel::load_default_prime()
{
    return StellarPopulationModel(IsochroneGrid::load_default_prime());
}

StellarPopulationModel::StellarPopulationModel(IsochroneGrid grid)
    : isochrones_(std::move(grid))
{
}

StellarProperties StellarPopulationModel::lookup(const StellarPopulationQuery &query) const
{
    std::string component = query.component;
    if (component.empty()) {
        component = component_name(query.component_index);
    }

    IsochroneQuery isochrone_query;
    isochrone_query.component = component;
    isochrone_query.initial_mass_msun = query.initial_mass_msun;
    isochrone_query.log_age = query.use_default_log_age ? default_log_age(component) : query.log_age;
    isochrone_query.metallicity_mh =
        query.use_default_metallicity ? default_metallicity_mh(component) : query.metallicity_mh;
    return isochrones_.lookup(isochrone_query);
}

std::vector<MassInterval> StellarPopulationModel::matching_initial_mass_intervals(
    const StellarPopulationQuery &query,
    const std::vector<MagnitudeSelection> &selection) const
{
    std::string component = query.component;
    if (component.empty()) {
        component = component_name(query.component_index);
    }

    IsochroneQuery isochrone_query;
    isochrone_query.component = component;
    isochrone_query.initial_mass_msun = query.initial_mass_msun;
    isochrone_query.log_age = query.use_default_log_age ? default_log_age(component) : query.log_age;
    isochrone_query.metallicity_mh =
        query.use_default_metallicity ? default_metallicity_mh(component) : query.metallicity_mh;
    return isochrones_.matching_initial_mass_intervals(isochrone_query, selection);
}

std::string StellarPopulationModel::component_name(int component_index)
{
    return find_component(component_index).name;
}

int StellarPopulationModel::component_index(const std::string &component)
{
    return find_component(component).index;
}

double StellarPopulationModel::default_log_age(const std::string &component)
{
    return find_component(component).log_age;
}

double StellarPopulationModel::default_metallicity_mh(const std::string &component)
{
    return find_component(component).metallicity_mh;
}

} // namespace genulens::model
