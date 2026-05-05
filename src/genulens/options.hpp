#pragma once

#include "genulens/model/parameters.hpp"
#include "genulens/simulation/observation_config.hpp"

#include <string>
#include <vector>

namespace genulens {

struct SourceSelectionConfig {
    double i_min = 14.0;
    double i_max = 21.0;
    double vi_min = 0.0;
    double vi_max = 0.0;
    double ai_rc = 0.0;
    double evi_rc = 0.0;
    double dm_rc = 0.0;
    double ak_rc = 0.0;
    double dust_scale_height_pc = 164.0;
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
    SamplingConfig sampling;
    RuntimeConfig runtime;
    std::string input_dir;
    std::vector<std::string> raw_cli_args;
};

} // namespace genulens
