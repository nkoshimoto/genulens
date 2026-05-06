#include "genulens/types.hpp"

#include <numeric>

namespace genulens {

double ComponentDensities::total() const
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

std::vector<std::string> SimulationResult::columns() const
{
    std::vector<std::string> out;
    if (verbosity == 1) {
        out = {"wtj", "tE", "thetaE", "piE", "M_L", "D_S", "D_L",
               "mu_rel", "iS", "iL", "tau_s", "tau_l", "fREM"};
    } else if (verbosity == 2) {
        out = {"wtj", "tE", "thetaE", "piEN", "piEE", "D_S",
               "muSl", "muSb", "iS", "iL", "fREM"};
    } else if (verbosity == 3 || verbosity == 0) {
        out = {"wtj", "M_L", "D_L", "D_S", "t_E", "theta_E",
               "pi_E", "pi_EN", "pi_EE", "mu_rel", "mu_rel_N", "mu_rel_E",
               "mu_Sl", "mu_Sb", "I_L", "K_L", "iS", "iL", "fREM"};
    } else if (verbosity == 4) {
        out = {"wtj", "M_L", "D_L", "D_S", "t_E", "theta_E",
               "pi_E", "pi_EN", "pi_EE", "mu_rel", "mu_rel_N", "mu_rel_E",
               "mu_Sl", "mu_Sb", "I_L", "K_L", "iS", "iL", "fREM"};
    } else if (verbosity >= 5) {
        out = {"wtj", "M_L", "D_L", "D_S", "t_E", "theta_E",
               "pi_E", "pi_EN", "pi_EE", "mu_rel", "mu_rel_N", "mu_rel_E",
               "mu_Sl", "mu_Sb", "I_L", "K_L", "iS", "iL", "fREM"};
    }
    if (include_source_properties) {
        out.push_back("logage_S");
        out.push_back("MH_S");
        out.push_back("Zini_S");
        out.push_back("M_S_ini");
        out.push_back("M_S");
        out.push_back("R_S");
        out.push_back("teff_S");
        out.push_back("logg_S");
        out.push_back("theta_S");
        for (const auto &band : source_property_bands) {
            out.push_back("M_" + band + "_S");
        }
    }
    return out;
}

std::vector<double> SimulationResult::flattened_rows() const
{
    std::vector<double> rows;
    rows.reserve(events.size() * column_count());
    for (const auto &event : events) {
        if (verbosity == 1) {
            rows.push_back(event.weight);
            rows.push_back(event.tE);
            rows.push_back(event.thetaE);
            rows.push_back(event.piE);
            rows.push_back(event.lens_mass_msun);
            rows.push_back(event.source_distance_pc);
            rows.push_back(event.lens_distance_pc);
            rows.push_back(event.mu_rel_masyr);
            rows.push_back(static_cast<double>(event.source_component));
            rows.push_back(static_cast<double>(event.lens_component));
            rows.push_back(event.source_age_gyr);
            rows.push_back(event.lens_age_gyr);
            rows.push_back(static_cast<double>(event.remnant_flag));
        } else if (verbosity == 2) {
            rows.push_back(event.weight);
            rows.push_back(event.tE);
            rows.push_back(event.thetaE);
            rows.push_back(event.piEN);
            rows.push_back(event.piEE);
            rows.push_back(event.source_distance_pc);
            rows.push_back(event.source_mu_l_masyr);
            rows.push_back(event.source_mu_b_masyr);
            rows.push_back(static_cast<double>(event.source_component));
            rows.push_back(static_cast<double>(event.lens_component));
            rows.push_back(static_cast<double>(event.remnant_flag));
        } else if (verbosity >= 3 || verbosity == 0) {
            rows.push_back(event.weight);
            rows.push_back(event.lens_mass_msun);
            rows.push_back(event.lens_distance_pc);
            rows.push_back(event.source_distance_pc);
            rows.push_back(event.tE);
            rows.push_back(event.thetaE);
            rows.push_back(event.piE);
            rows.push_back(event.piEN);
            rows.push_back(event.piEE);
            rows.push_back(event.mu_rel_masyr);
            rows.push_back(event.mu_rel_N_masyr);
            rows.push_back(event.mu_rel_E_masyr);
            rows.push_back(event.source_mu_l_masyr);
            rows.push_back(event.source_mu_b_masyr);
            rows.push_back(event.lens_i_mag);
            rows.push_back(event.lens_k_mag);
            rows.push_back(static_cast<double>(event.source_component));
            rows.push_back(static_cast<double>(event.lens_component));
            rows.push_back(static_cast<double>(event.remnant_flag));
        }
        if (include_source_properties) {
            rows.push_back(event.source_log_age);
            rows.push_back(event.source_metallicity_mh);
            rows.push_back(event.source_zini);
            rows.push_back(event.source_initial_mass_msun);
            rows.push_back(event.source_current_mass_msun);
            rows.push_back(event.source_radius_rsun);
            rows.push_back(event.source_teff_k);
            rows.push_back(event.source_logg);
            rows.push_back(event.source_angular_radius_microarcsec);
            for (std::size_t i = 0; i < source_property_bands.size(); ++i) {
                rows.push_back(i < event.source_absolute_magnitudes.size()
                                   ? event.source_absolute_magnitudes[i]
                                   : Event::missing_value());
            }
        }
    }
    return rows;
}

} // namespace genulens
