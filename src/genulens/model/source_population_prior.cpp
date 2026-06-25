#include "genulens/model/source_population_prior.hpp"

#include "genulens/model/stellar_population_model.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace genulens::model {
namespace {

AgeMetallicityPoint point(double log_age, double metallicity_mh, double weight)
{
    return {log_age, metallicity_mh, weight};
}

struct WeightedAge {
    double log_age = 0.0;
    double weight = 1.0;
};

struct WeightedMetallicity {
    double metallicity_mh = 0.0;
    double weight = 1.0;
};

double log_age_from_gyr(double age_gyr)
{
    return 9.0 + std::log10(age_gyr);
}

WeightedAge age_gyr(double age_gyr, double weight)
{
    return {log_age_from_gyr(age_gyr), weight};
}

std::vector<AgeMetallicityPoint> combine(
    const std::vector<WeightedAge> &ages,
    const std::vector<WeightedMetallicity> &metallicities)
{
    double total_weight = 0.0;
    for (const auto &age : ages) {
        for (const auto &metallicity : metallicities) {
            total_weight += age.weight * metallicity.weight;
        }
    }
    if (total_weight <= 0.0) {
        throw std::runtime_error("source-population prior has non-positive total weight");
    }

    std::vector<AgeMetallicityPoint> points;
    points.reserve(ages.size() * metallicities.size());
    for (const auto &age : ages) {
        for (const auto &metallicity : metallicities) {
            points.push_back(point(age.log_age,
                                   metallicity.metallicity_mh,
                                   age.weight * metallicity.weight / total_weight));
        }
    }
    return points;
}

std::vector<WeightedAge> disk_age_mixture(int component_index)
{
    // These are a compact quadrature over the age ranges used by the legacy
    // genstars/Koshimoto luminosity-function tables, in Gyr.
    static constexpr std::array<std::array<double, 3>, 7> kAgeGyr = {{
        {{0.05, 0.10, 0.15}},
        {{0.15, 0.586449, 1.0}},
        {{1.0, 1.516357, 2.0}},
        {{2.0, 2.516884, 3.0}},
        {{3.0, 4.068387, 5.0}},
        {{5.0, 6.069263, 7.0}},
        {{7.0, 8.656024, 10.0}},
    }};
    static constexpr double kThinDiskSfrTimescaleGyr = 7.0;

    std::vector<WeightedAge> ages;
    ages.reserve(3);
    for (std::size_t i = 0; i < 3; ++i) {
        const double age = kAgeGyr.at(component_index).at(i);
        double weight = std::exp(-(10.0 - age) / kThinDiskSfrTimescaleGyr);
        if (i == 0 || i == 2) weight *= 0.5;
        ages.push_back(age_gyr(age, weight));
    }
    return ages;
}

std::vector<WeightedAge> gaussian_age_mixture(double center_gyr, double sigma_gyr)
{
    const std::array<double, 3> ages = {
        center_gyr - sigma_gyr,
        center_gyr,
        center_gyr + sigma_gyr,
    };

    std::vector<WeightedAge> points;
    points.reserve(ages.size());
    for (const double age : ages) {
        const double z = (age - center_gyr) / sigma_gyr;
        points.push_back(age_gyr(age, std::exp(-0.5 * z * z)));
    }
    return points;
}

std::vector<WeightedMetallicity> disk_metallicity_mixture(int component_index)
{
    if (component_index <= 2) {
        return {
            {-0.5, 0.10},
            {-0.25, 0.20},
            {0.0, 0.40},
            {0.25, 0.30},
        };
    }
    if (component_index <= 4) {
        return {
            {-0.5, 0.15},
            {-0.25, 0.30},
            {0.0, 0.35},
            {0.25, 0.20},
        };
    }
    return {
        {-0.5, 0.25},
        {-0.25, 0.35},
        {0.0, 0.30},
        {0.25, 0.10},
    };
}

} // namespace

std::vector<AgeMetallicityPoint> SourcePopulationPrior::points_for_component(int component_index)
{
    if (component_index >= 0 && component_index <= 6) {
        return combine(disk_age_mixture(component_index),
                       disk_metallicity_mixture(component_index));
    }
    if (component_index == 7 || component_index == 10) {
        return combine({age_gyr(12.0, 1.0)},
                       {{-1.2, 0.10}, {-1.0, 0.20}, {-0.8, 0.40}, {-0.6, 0.20}, {-0.4, 0.10}});
    }
    if (component_index == 8) {
        return combine(gaussian_age_mixture(9.0, 1.0),
                       {{-0.5, 0.20}, {-0.25, 0.25}, {0.0, 0.30}, {0.25, 0.25}});
    }
    if (component_index == 9) {
        return combine(gaussian_age_mixture(7.0, 1.0),
                       {{-0.5, 0.10}, {-0.25, 0.25}, {0.0, 0.40}, {0.25, 0.25}});
    }
    throw std::runtime_error("unknown source-population prior component index: " +
                             std::to_string(component_index));
}

} // namespace genulens::model
