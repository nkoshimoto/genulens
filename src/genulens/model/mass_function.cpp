#include "genulens/model/mass_function.hpp"

#include "genulens/model/parameters.hpp"

#include <cmath>

namespace genulens::model {

namespace {

constexpr double kNeutronStarMinMass = 1.2;
constexpr double kNeutronStarMaxMass = 2.1;
constexpr double kMaxSigmaLogA = 1.8;
constexpr double kMinSigmaLogA = 0.3;
constexpr double kMaxMeanLogA = 1.7;
constexpr double kMinMeanLogA = 0.6;

} // namespace

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
    const double temp44 = std::pow(p.mbr, p.alpha4 + 1.0);
    const double temp45 = std::pow(p.mbr, p.alpha5 + 1.0);

    for (int i = 0; i <= n_bins; ++i) {
        const double log_mass = i * grid.log_mass_step + grid.log_mass_start;
        grid.log_mass[i] = log_mass;
        grid.mass[i] = std::pow(10.0, log_mass);
        const double mass = grid.mass[i];
        const double alpha = (mass < p.mbr) ? p.alpha5 :
                             (mass < p.m3) ? p.alpha4 :
                             (mass < p.m2) ? p.alpha3 :
                             (mass < p.m1) ? p.alpha2 :
                             (mass < p.m0) ? p.alpha1 : p.alpha0;
        double continuity = 1.0;
        if (mass < p.m0) continuity = temp01 / temp00;
        if (mass < p.m1) continuity = temp12 / temp11 * continuity;
        if (mass < p.m2) continuity = temp23 / temp22 * continuity;
        if (mass < p.m3) continuity = temp34 / temp33 * continuity;
        if (mass < p.mbr) continuity = temp45 / temp44 * continuity;
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

RemnantMassModel::RemnantMassModel(double white_dwarf_initial_mass_max)
    : white_dwarf_initial_mass_max_(white_dwarf_initial_mass_max)
{
}

RemnantMass RemnantMassModel::evolve(double initial_mass_msun, bool mean, genulens::RandomEngine &rng) const
{
    RemnantMass remnant;
    const double mini = initial_mass_msun;
    const double p_ns = (mini < white_dwarf_initial_mass_max_) ? 0.0 :
                        (mini < 15.0) ? 1.0 :
                        (mini < 17.8) ? 0.679 :
                        (mini < 18.5) ? 0.833 :
                        (mini < 21.7) ? 0.500 :
                        (mini < 25.2) ? 0.0 :
                        (mini < 27.5) ? 0.652 :
                        (mini < 60.0) ? 0.0 : 0.4;
    if (mini < white_dwarf_initial_mass_max_) {
        remnant.mass_msun = 0.109 * mini + 0.394;
        remnant.remnant_type = 1.0;
        return remnant;
    }

    double m_ns = 0.0;
    do {
        m_ns = (mini < 13.0) ? 2.24 + 0.508 * (mini - 14.75)
                                     + 0.125 * (mini - 14.75) * (mini - 14.75)
                                     + 0.011 * (mini - 14.75) * (mini - 14.75) * (mini - 14.75)
             : (mini < 15.0) ? 0.123 + 0.112 * mini
             : (mini < 17.8) ? 0.996 + 0.0384 * mini
             : (mini < 18.5) ? -0.020 + 0.10 * mini
             : (mini < 21.7 && !mean) ? 1.60 + 0.158 * rng.gaussian()
             : (mini < 21.7 && mean) ? 1.60
             : (mini < 27.5) ? 3232.29 - 409.429 * (mini - 2.619)
                                       + 17.2867 * (mini - 2.619) * (mini - 2.619)
                                       - 0.24315 * (mini - 2.619) * (mini - 2.619) * (mini - 2.619)
             : (!mean) ? 1.78 + 0.02 * rng.gaussian()
             : 1.78;
    } while (p_ns > 0.0 && (m_ns < kNeutronStarMinMass || m_ns > kNeutronStarMaxMass));

    const double m_core = (mini < 42.21) ? -2.049 + 0.4140 * mini
                                         : 5.697 + 7.8598e8 * std::pow(mini, -4.858);
    const double m_all = 15.52 - 0.3294 * (mini - 25.97)
                         - 0.02121 * (mini - 25.97) * (mini - 25.97)
                         + 0.003120 * (mini - 25.97) * (mini - 25.97) * (mini - 25.97);
    const double f_ej = (mini < 42.21) ? 0.9 : 1.0;
    const double m_bh = f_ej * m_core + (1.0 - f_ej) * m_all;

    if (mean) {
        remnant.mass_msun = p_ns * m_ns + (1.0 - p_ns) * m_bh;
        remnant.remnant_type = p_ns * 2.0 + (1.0 - p_ns) * 3.0;
    } else {
        const double draw = rng.uniform();
        remnant.mass_msun = (draw < p_ns) ? m_ns : m_bh;
        remnant.remnant_type = (draw < p_ns) ? 2.0 : 3.0;
    }
    return remnant;
}

ProjectedSeparation BinaryLensSampler::sample_projected_separation(double primary_mass_msun,
                                                                    double secondary_mass_msun,
                                                                    int coefficient,
                                                                    genulens::RandomEngine &rng) const
{
    const double primary = (secondary_mass_msun > primary_mass_msun) ? secondary_mass_msun : primary_mass_msun;
    double mean_log_a = 0.57 + 1.02 * primary;
    if (mean_log_a > kMaxMeanLogA) mean_log_a = kMaxMeanLogA;
    if (mean_log_a < kMinMeanLogA) mean_log_a = kMinMeanLogA;

    double sigma_log_a = 1.61 + 1.15 * std::log10(primary);
    if (sigma_log_a > kMaxSigmaLogA) sigma_log_a = kMaxSigmaLogA;
    if (sigma_log_a < kMinSigmaLogA) sigma_log_a = kMinSigmaLogA;

    const double signed_half_normal = coefficient * std::fabs(rng.gaussian());
    const double log_a = mean_log_a + signed_half_normal * sigma_log_a;
    const double semi_major_axis = std::pow(10.0, log_a);
    const double projection_draw = rng.uniform();

    ProjectedSeparation result;
    result.log_separation_au = log_a;
    result.projected_separation_au = std::sqrt(1.0 - projection_draw * projection_draw) * semi_major_axis;
    return result;
}

double broken_power_law_imf(double mass_msun, const IMFParameters &params)
{
    if (mass_msun <= 0.0) return 0.0;
    if (mass_msun < params.mbr) return std::pow(mass_msun, params.alpha5);
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
