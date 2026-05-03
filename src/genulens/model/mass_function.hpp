#pragma once

#include "genulens/model/parameters.hpp"

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

class BrokenPowerLawIMF {
public:
    explicit BrokenPowerLawIMF(IMFParameters parameters);

    double evaluate(double mass_msun) const;
    MassFunctionGrid build_grid(int n_bins, double lower_mass_msun, double upper_mass_msun) const;

private:
    IMFParameters parameters_;
};

double broken_power_law_imf(double mass_msun, const IMFParameters &params);
double broken_power_law_imf(double mass_msun);

} // namespace genulens::model
