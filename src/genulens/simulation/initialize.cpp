#include "genulens/simulation/initialize.hpp"

#include "genulens/cli/option.h"

namespace genulens {

RunContext Initializer::create_context() const
{
    return {};
}

void Initializer::read_sampling_options(RunContext &context, int argc, char **argv,
                                                double cos_pa, double sin_pa) const
{
    auto &options = context.sampling;
    options.n_simu = getOptiond(argc, argv, "Nsimu", 1, options.n_simu);
    options.n_like_min = getOptiond(argc, argv, "NlikeMIN", 1, options.n_like_min);
    options.v_earth_l = getOptiond(argc, argv, "vEarthlb", 1, options.v_earth_l);
    options.v_earth_b = getOptiond(argc, argv, "vEarthlb", 2, options.v_earth_b);
    const double v_earth_e = getOptiond(argc, argv, "vEarthEN", 1, 0);
    const double v_earth_n = getOptiond(argc, argv, "vEarthEN", 2, 0);
    if (options.v_earth_l == 11.9392 && options.v_earth_b == -17.9020 &&
        v_earth_n != 0 && v_earth_e != 0) {
        options.v_earth_b = cos_pa * v_earth_n - sin_pa * v_earth_e;
        options.v_earth_l = sin_pa * v_earth_n + cos_pa * v_earth_e;
    }
    options.gamma_ds = getOptiond(argc, argv, "gammaDs", 1, options.gamma_ds);
    options.weight_lens_distance = getOptiond(argc, argv, "wtD_L", 1, options.weight_lens_distance);
    options.weight_lens_mass = getOptiond(argc, argv, "wtM_L", 1, options.weight_lens_mass);
    options.no_gamma_importance_sampling =
        getOptiond(argc, argv, "NoGAMMAIS", 1, options.no_gamma_importance_sampling);
    options.small_gamma = getOptiond(argc, argv, "SMALLGAMMA", 1, options.small_gamma);
    options.verbosity = getOptiond(argc, argv, "VERBOSITY", 1, options.verbosity);
    options.uniform_likelihood = getOptiond(argc, argv, "UNIFORM", 1, options.uniform_likelihood);
    options.binary = getOptiond(argc, argv, "BINARY", 1, options.binary);
    options.remnant = getOptiond(argc, argv, "REMNANT", 1, options.remnant);
    options.only_white_dwarf = getOptiond(argc, argv, "onlyWD", 1, options.only_white_dwarf);
    if (options.only_white_dwarf == 1) options.remnant = 0;
    options.calc_prior_piE = getOptiond(argc, argv, "CALCPRIORpiE", 1, options.calc_prior_piE);
    options.calc_prior_thetaE = getOptiond(argc, argv, "CALCPRIORthE", 1, options.calc_prior_thetaE);
}

} // namespace genulens
