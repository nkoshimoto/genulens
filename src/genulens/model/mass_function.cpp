#include "genulens/model/mass_function.hpp"

#include "genulens/model/parameters.hpp"

#include <cmath>

namespace genulens::model {

BrokenPowerLawIMF::BrokenPowerLawIMF(IMFParameters parameters)
    : parameters_(parameters)
{
}

double BrokenPowerLawIMF::evaluate(double mass_msun) const
{
    return broken_power_law_imf(mass_msun, parameters_);
}

MassFunctionGrid BrokenPowerLawIMF::build_grid(int n_bins, double lower_mass_msun, double upper_mass_msun) const
{
    MassFunctionGrid grid;
    grid.log_mass.resize(n_bins + 1);
    grid.mass.resize(n_bins + 1);
    grid.probability_log_mass.resize(n_bins + 1);
    grid.cumulative_number.resize(n_bins + 1);
    grid.cumulative_number_norm.resize(n_bins + 1);
    grid.cumulative_mass.resize(n_bins + 1);
    grid.cumulative_mass_norm.resize(n_bins + 1);
    grid.percentile_index.assign(22, 0);
    grid.log_mass_start = std::log10(lower_mass_msun);
    grid.log_mass_step = (std::log10(upper_mass_msun) - grid.log_mass_start) / n_bins;

    const auto &p = parameters_;
    const double temp00 = std::pow(p.m0, p.alpha0 + 1.0);
    const double temp01 = std::pow(p.m0, p.alpha1 + 1.0);
    const double temp11 = std::pow(p.m1, p.alpha1 + 1.0);
    const double temp12 = std::pow(p.m1, p.alpha2 + 1.0);
    const double temp22 = std::pow(p.m2, p.alpha2 + 1.0);
    const double temp23 = std::pow(p.m2, p.alpha3 + 1.0);
    const double temp33 = std::pow(p.m3, p.alpha3 + 1.0);
    const double temp34 = std::pow(p.m3, p.alpha4 + 1.0);

    for (int i = 0; i <= n_bins; ++i) {
        const double log_mass = i * grid.log_mass_step + grid.log_mass_start;
        grid.log_mass[i] = log_mass;
        grid.mass[i] = std::pow(10.0, log_mass);
        const double mass = grid.mass[i];
        const double alpha = (mass < p.m3) ? p.alpha4 :
                             (mass < p.m2) ? p.alpha3 :
                             (mass < p.m1) ? p.alpha2 :
                             (mass < p.m0) ? p.alpha1 : p.alpha0;
        double continuity = 1.0;
        if (mass < p.m0) continuity = temp01 / temp00;
        if (mass < p.m1) continuity = temp12 / temp11 * continuity;
        if (mass < p.m2) continuity = temp23 / temp22 * continuity;
        if (mass < p.m3) continuity = temp34 / temp33 * continuity;
        grid.probability_log_mass[i] = std::pow(mass, alpha + 1.0) / continuity;
        if (i >= 1) {
            grid.cumulative_number[i] =
                0.5 * (grid.probability_log_mass[i] + grid.probability_log_mass[i - 1]) *
                grid.log_mass_step + grid.cumulative_number[i - 1];
            grid.cumulative_mass[i] =
                0.5 * (grid.mass[i] * grid.probability_log_mass[i] +
                       grid.mass[i - 1] * grid.probability_log_mass[i - 1]) *
                grid.log_mass_step + grid.cumulative_mass[i - 1];
        }
    }

    for (int i = 0; i <= n_bins; ++i) {
        grid.cumulative_number_norm[i] = grid.cumulative_number[i] / grid.cumulative_number[n_bins];
        grid.cumulative_mass_norm[i] = grid.cumulative_mass[i] / grid.cumulative_mass[n_bins];
        grid.probability_log_mass[i] /= grid.cumulative_number[n_bins];
        const int percentile = static_cast<int>(grid.cumulative_number_norm[i] * 20);
        if (grid.percentile_index[percentile] == 0) {
            grid.percentile_index[percentile] = (percentile == 0) ? 1 : static_cast<int>(i + 0.5);
        }
    }
    return grid;
}

double broken_power_law_imf(double mass_msun, const IMFParameters &params)
{
    if (mass_msun <= 0.0) return 0.0;
    if (mass_msun < params.m3) return std::pow(mass_msun, params.alpha4);
    if (mass_msun < params.m2) return std::pow(mass_msun, params.alpha3);
    if (mass_msun < params.m1) return std::pow(mass_msun, params.alpha2);
    if (mass_msun < params.m0) return std::pow(mass_msun, params.alpha1);
    return std::pow(mass_msun, params.alpha0);
}

double broken_power_law_imf(double mass_msun)
{
    return broken_power_law_imf(mass_msun, default_model_parameters().imf);
}

} // namespace genulens::model
