#pragma once

#include "genulens/model/parameters.hpp"
#include "genulens/rng.hpp"

#include <vector>

namespace genulens::model {

struct MassFunctionGrid {
    std::vector<double> log_mass;
    std::vector<double> mass;
    std::vector<double> probability_log_mass;
    std::vector<double> cumulative_number;
    std::vector<double> cumulative_number_norm;
    std::vector<double> cumulative_mass;
    std::vector<double> cumulative_mass_norm;
    std::vector<int> percentile_index;
    double log_mass_start = 0.0;
    double log_mass_step = 0.0;
};

struct RemnantMass {
    double mass_msun = 0.0;
    double remnant_type = 0.0;
};

struct ProjectedSeparation {
    double log_separation_au = 0.0;
    double projected_separation_au = 0.0;
};

class BrokenPowerLawIMF {
public:
    explicit BrokenPowerLawIMF(IMFParameters parameters);

    double evaluate(double mass_msun) const;
    MassFunctionGrid build_grid(int n_bins, double lower_mass_msun, double upper_mass_msun) const;

private:
    IMFParameters parameters_;
};

class RemnantMassModel {
public:
    explicit RemnantMassModel(double white_dwarf_initial_mass_max = 9.0);

    RemnantMass evolve(double initial_mass_msun, bool mean, genulens::RandomEngine &rng) const;

private:
    double white_dwarf_initial_mass_max_ = 9.0;
};

class BinaryLensSampler {
public:
    ProjectedSeparation sample_projected_separation(double primary_mass_msun,
                                                     double secondary_mass_msun,
                                                     int coefficient,
                                                     genulens::RandomEngine &rng) const;
};

double broken_power_law_imf(double mass_msun, const IMFParameters &params);
double broken_power_law_imf(double mass_msun);

} // namespace genulens::model
