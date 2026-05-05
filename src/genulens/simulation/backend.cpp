#include "genulens/simulation/backend.hpp"

#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/sampler.hpp"
#include "genulens/rng.hpp"

#include <memory>
#include <utility>

namespace genulens {

namespace {

bool is_typed_config_key(const std::string &arg)
{
    return arg == "VERBOSITY" || arg == "l" || arg == "b" || arg == "Nsimu" || arg == "seed" ||
           arg == "tE" || arg == "M0" || arg == "M1" || arg == "M2" || arg == "M3" ||
           arg == "Ml" || arg == "Mu" || arg == "alpha0" || arg == "alpha1" ||
           arg == "alpha2" || arg == "alpha3" || arg == "alpha4" ||
           arg == "Isrange" || arg == "VIsrange" || arg == "AIrc" || arg == "EVIrc" ||
           arg == "DMrc" || arg == "AKrc" || arg == "hdust" ||
           arg == "NlikeMIN" || arg == "vEarthlb" || arg == "gammaDs" ||
           arg == "wtD_L" || arg == "wtM_L" || arg == "NoGAMMAIS" ||
           arg == "SMALLGAMMA" || arg == "UNIFORM" || arg == "BINARY" ||
           arg == "REMNANT" || arg == "onlyWD" || arg == "CALCPRIORpiE" ||
           arg == "CALCPRIORthE" || arg == "Dmax" || arg == "CALCTAU" ||
           arg == "thetaE" || arg == "thetaEdet" || arg == "thetaErange" ||
           arg == "piE" || arg == "piEdet" || arg == "piErange" ||
           arg == "piEN" || arg == "piEE" || arg == "musl" || arg == "musb" ||
           arg == "musN" || arg == "musE" || arg == "musRCG" ||
           arg == "muhelN" || arg == "muhelE" || arg == "IL" || arg == "ILdet" ||
           arg == "KL" || arg == "KLdet" || arg == "u0" ||
           arg == "DISK" || arg == "rhot0" || arg == "hDISK" || arg == "addX" ||
           arg == "model" || arg == "R0" || arg == "thetaD" || arg == "frho0b" ||
           arg == "Rc" || arg == "zb_c" || arg == "x0" || arg == "y0" ||
           arg == "z0" || arg == "C1" || arg == "C2" || arg == "C3" ||
           arg == "x0_X" || arg == "y0_X" || arg == "z0_X" || arg == "C1_X" ||
           arg == "C2_X" || arg == "b_zX" || arg == "fX" || arg == "Rc_X" ||
           arg == "b_zY" || arg == "SH" || arg == "rho0SHMS" ||
           arg == "Omega_p" || arg == "model_vb" || arg == "x0_vb" ||
           arg == "y0_vb" || arg == "z0_vb" || arg == "C1_vb" ||
           arg == "C2_vb" || arg == "C3_vb" || arg == "sigx_vb" ||
           arg == "sigy_vb" || arg == "sigz_vb" || arg == "sigx_vb0" ||
           arg == "sigy_vb0" || arg == "sigz_vb0" || arg == "vx_str" ||
           arg == "y0_str" || arg == "model_vbz" || arg == "x0_vbz" ||
           arg == "y0_vbz" || arg == "z0_vbz" || arg == "C1_vbz" ||
           arg == "C2_vbz" || arg == "C3_vbz" || arg == "hsigUt" ||
           arg == "hsigWt" || arg == "hsigUT" || arg == "hsigWT" ||
           arg == "betaU" || arg == "betaW" || arg == "sigU10d" ||
           arg == "sigW10d" || arg == "sigU0td" || arg == "sigW0td" ||
           arg == "sigU_SH" || arg == "sigV_SH" || arg == "sigW_SH" ||
           arg == "NSD" || arg == "x0ND" || arg == "y0ND" ||
           arg == "z0ND" || arg == "MND" || arg == "MXDkick" ||
           arg == "vkickBH" || arg == "vkickNS" || arg == "BHhd" ||
           arg == "BHhb" || arg == "fixRhdBH" || arg == "RhdBH0" ||
           arg == "betaBH" || arg == "UseSigBH";
}

int typed_config_value_count(const std::string &arg)
{
    if (arg == "tE" || arg == "Isrange" || arg == "VIsrange" || arg == "vEarthlb" ||
        arg == "thetaE" || arg == "thetaErange" || arg == "piE" || arg == "piErange" ||
        arg == "piEN" || arg == "piEE" || arg == "musl" || arg == "musb" ||
        arg == "musN" || arg == "musE" || arg == "muhelN" || arg == "muhelE" ||
        arg == "IL" || arg == "KL") {
        return 2;
    }
    return 1;
}

std::vector<std::string> legacy_args_without_typed_keys(const GenulensConfig &config)
{
    std::vector<std::string> args;
    args.push_back(config.raw_cli_args.empty() ? "genulens" : config.raw_cli_args.front());
    for (std::size_t i = 1; i < config.raw_cli_args.size(); ++i) {
        const auto &arg = config.raw_cli_args[i];
        if (is_typed_config_key(arg)) {
            i += typed_config_value_count(arg);
            continue;
        }
        args.push_back(arg);
    }
    return args;
}

} // namespace

SimulationResult SimulationBackend::simulate(const GenulensConfig &config, LikelihoodFunction likelihood) const
{
    auto args = legacy_args_without_typed_keys(config);
    std::vector<char *> argv;
    argv.reserve(args.size());
    for (const auto &arg : args) {
        argv.push_back(const_cast<char *>(arg.c_str()));
    }

    Initializer initializer;
    auto context = initializer.create_context();
    context.seed = config.seed;
    context.runtime.rng = std::make_unique<RandomEngine>(config.seed);
    return Sampler().simulate(context, config, static_cast<int>(argv.size()), argv.data(), std::move(likelihood));
}

} // namespace genulens
