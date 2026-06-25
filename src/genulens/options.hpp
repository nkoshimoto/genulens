#pragma once

#include "genulens/model/parameters.hpp"
#include "genulens/simulation/observation_config.hpp"

#include <string>
#include <vector>

namespace genulens {

struct SourceSelectionConfig {
    std::string mode = "classic";
    std::string photometry = "prime";
    std::string band = "Imag";
    double min_magnitude = 14.0;
    double max_magnitude = 21.0;
    int apparent_magnitude = 1;
    int use_magnitude_selection = 1;

    double i_min = 14.0;
    double i_max = 21.0;
    double vi_min = 0.0;
    double vi_max = 0.0;
    double av_rc = 0.0;
    double ai_rc = 0.0;
    double aj_rc = 0.0;
    double ah_rc = 0.0;
    double evi_rc = 0.0;
    double dm_rc = 0.0;
    double ak_rc = 0.0;
    std::string extinction_mode = "manual";
    double ejk_rc = 0.0;
    double ejk_scale = 1.0;
    double ejk_offset = 0.0;
    int extinction_law = 1;
    int extinction_map = 1;
    std::string extinction_map_path;
    double dust_scale_height_pc = 164.0;

    double min_initial_mass_msun = 0.09;
    double max_initial_mass_msun = 1.0;
    std::string isochrone_model = "parsec_solar_scaled";
    std::string isochrone_family = "parsec";
    std::string isochrone_abundance = "solar_scaled";
    double isochrone_alpha_fe = 0.0;
    std::string isochrone_table_path;
    std::string secondary_isochrone_family = "mist";
    std::string secondary_isochrone_abundance = "alpha_enhanced";
    double secondary_isochrone_alpha_fe = 0.4;
    std::string secondary_isochrone_table_path;
    double alpha_enhanced_fraction = 0.0;
    std::vector<int> alpha_enhanced_components;
    std::vector<double> alpha_enhanced_component_fractions;
};

struct ForwardSourceConfig {
    int enabled = 0;
    int use_model_imf = 1;
    model::IMFParameters imf;
    std::string photometry = "roman";
    std::string isochrone_model = "parsec_solar_scaled";
    std::string isochrone_family = "parsec";
    std::string isochrone_abundance = "solar_scaled";
    double isochrone_alpha_fe = 0.0;
    std::string isochrone_table_path;
    std::string secondary_isochrone_family = "mist";
    std::string secondary_isochrone_abundance = "alpha_enhanced";
    double secondary_isochrone_alpha_fe = 0.4;
    std::string secondary_isochrone_table_path;
    std::string alpha_enhanced_table_path;
    double alpha_enhanced_fraction = 0.0;
    std::vector<int> alpha_enhanced_components;
    std::vector<double> alpha_enhanced_component_fractions;
    double min_initial_mass_msun = 0.09;
    double max_initial_mass_msun = 1.0;
    std::vector<std::string> selection_bands;
    std::vector<double> selection_min_magnitudes;
    std::vector<double> selection_max_magnitudes;
    std::vector<int> selection_apparent_magnitudes;
};

struct SamplingConfig {
    long n_like_min = 0;
    double v_earth_l = 11.9392;
    double v_earth_b = -17.9020;
    double gamma_ds = 0.5;
    double weight_lens_distance = 0.0;
    double weight_lens_mass = 0.0;
    int no_gamma_importance_sampling = 0;
    int small_gamma = 0;
    int verbosity = 0;
    int uniform_likelihood = 0;
    int binary = 0;
    int remnant = 0;
    int only_white_dwarf = 0;
    int calc_prior_piE = 0;
    int calc_prior_thetaE = 0;
};

struct RuntimeConfig {
    int max_distance_pc = 16000;
    int calculate_optical_depth = 0;
};

struct GenulensConfig {
    unsigned long seed = 12304357UL;
    double l = 1.0;
    double b = -3.9;
    long n_simu = 100000;
    double observed_tE = 54.5;
    double observed_tE_error = 99999999999.0;
    bool use_gaussian_likelihood = true;
    model::ModelParameters model;
    ObservationConfig observation;
    SourceSelectionConfig source;
    ForwardSourceConfig forward_source;
    SamplingConfig sampling;
    RuntimeConfig runtime;
    std::string input_dir;
    std::vector<std::string> raw_cli_args;
};

} // namespace genulens
