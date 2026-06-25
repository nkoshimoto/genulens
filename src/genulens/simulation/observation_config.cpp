#include "genulens/simulation/observation_config.hpp"
#include "genulens/cli/option.h"
#include <cmath>

namespace genulens {

ObservationConfig ObservationConfig::from_cli(int argc, char **argv, int uniform_likelihood) {
    ObservationConfig obs;

    obs.tE_obs     = getOptiond(argc, argv, "tE",       1, 54.5);
    obs.tE_err     = getOptiond(argc, argv, "tE",       2, 99999999999.0);
    obs.tE_det     = (int)getOptiond(argc, argv, "tEdet",    1, 0);
    obs.tE_min     = getOptiond(argc, argv, "tErange",  1, 0);
    obs.tE_max     = getOptiond(argc, argv, "tErange",  2, 0);

    obs.thetaE_obs = getOptiond(argc, argv, "thetaE",     1, 0);
    obs.thetaE_err = getOptiond(argc, argv, "thetaE",     2, 0);
    obs.thetaE_det = (int)getOptiond(argc, argv, "thetaEdet",  1, 0);
    obs.thetaE_min = getOptiond(argc, argv, "thetaErange", 1, 0);
    obs.thetaE_max = getOptiond(argc, argv, "thetaErange", 2, 0);

    obs.piE_obs    = getOptiond(argc, argv, "piE",     1, 0);
    obs.piE_err    = getOptiond(argc, argv, "piE",     2, 0);
    obs.piE_det    = (int)getOptiond(argc, argv, "piEdet",  1, 0);
    obs.piE_min    = getOptiond(argc, argv, "piErange", 1, 0);
    obs.piE_max    = getOptiond(argc, argv, "piErange", 2, 0);

    obs.piEN_obs   = getOptiond(argc, argv, "piEN", 1, 0);
    obs.piEN_err   = getOptiond(argc, argv, "piEN", 2, 0);
    obs.piEE_obs   = getOptiond(argc, argv, "piEE", 1, 0);
    obs.piEE_err   = getOptiond(argc, argv, "piEE", 2, 0);

    obs.musl_obs   = getOptiond(argc, argv, "musl", 1, 0);
    obs.musl_err   = getOptiond(argc, argv, "musl", 2, 0);
    obs.musb_obs   = getOptiond(argc, argv, "musb", 1, 0);
    obs.musb_err   = getOptiond(argc, argv, "musb", 2, 0);
    obs.musN_obs   = getOptiond(argc, argv, "musN", 1, 0);
    obs.musN_err   = getOptiond(argc, argv, "musN", 2, 0);
    obs.musE_obs   = getOptiond(argc, argv, "musE", 1, 0);
    obs.musE_err   = getOptiond(argc, argv, "musE", 2, 0);
    obs.musRCG     = (int)getOptiond(argc, argv, "musRCG", 1, 0);

    obs.muhelN_obs = getOptiond(argc, argv, "muhelN", 1, 0);
    obs.muhelN_err = getOptiond(argc, argv, "muhelN", 2, 0);
    obs.muhelE_obs = getOptiond(argc, argv, "muhelE", 1, 0);
    obs.muhelE_err = getOptiond(argc, argv, "muhelE", 2, 0);

    obs.IL_obs     = getOptiond(argc, argv, "IL",    1, 14.00);
    obs.IL_err     = getOptiond(argc, argv, "IL",    2,  0.01);
    obs.IL_det     = (int)getOptiond(argc, argv, "ILdet", 1, 2);
    obs.KL_obs     = getOptiond(argc, argv, "KL",    1, 0);
    obs.KL_err     = getOptiond(argc, argv, "KL",    2, 0);
    obs.KL_det     = (int)getOptiond(argc, argv, "KLdet", 1, 0);
    obs.u0_obs     = getOptiond(argc, argv, "u0",    1, 0);

    // Importance sampling range setup
    int    NOIS      = (int)getOptiond(argc, argv, "NOIS", 1, 0);
    double fIS0      = (uniform_likelihood == 1) ? 1.02 : 4.0;
    double fIStE     = getOptiond(argc, argv, "fIStE",     1, fIS0);
    double fISthetaE = getOptiond(argc, argv, "fISthetaE", 1, fIS0);
    double fISpiE    = getOptiond(argc, argv, "fISpiE",    1, fIS0);

    if (NOIS == 0) {
        if (obs.tE_obs - obs.tE_err > 0 && obs.tE_err > 0 &&
            obs.tE_max - obs.tE_min == 0) {
            obs.tE_min = obs.tE_obs - fIStE * obs.tE_err;
            obs.tE_max = obs.tE_obs + fIStE * obs.tE_err;
            if (obs.tE_min <= 0 || obs.tE_det == 1) obs.tE_min = 1e-10;
            if (obs.tE_det == 2) obs.tE_max = 1e+6;
        }
        if (obs.thetaE_obs - obs.thetaE_err > 0 && obs.thetaE_err > 0 &&
            obs.thetaE_max - obs.thetaE_min == 0) {
            obs.thetaE_min = obs.thetaE_obs - fISthetaE * obs.thetaE_err;
            obs.thetaE_max = obs.thetaE_obs + fISthetaE * obs.thetaE_err;
            if (obs.thetaE_min <= 0 || obs.thetaE_det == 1) obs.thetaE_min = 1e-10;
            if (obs.thetaE_det == 2) obs.thetaE_max = 1e+6;
        }
        if (obs.piE_obs - obs.piE_err > 0 && obs.piE_err > 0 &&
            obs.piE_max - obs.piE_min == 0) {
            obs.piE_min = obs.piE_obs - fISpiE * obs.piE_err;
            obs.piE_max = obs.piE_obs + fISpiE * obs.piE_err;
            if (obs.piE_min <= 0 || obs.piE_det == 1) obs.piE_min = 1e-10;
            if (obs.piE_det == 2) obs.piE_max = 1e+6;
        }
        if (obs.piEN_err > 0 && obs.piEE_err > 0 &&
            obs.piE_max - obs.piE_min == 0) {
            double pNmin = fabs(obs.piEN_obs) - fISpiE * obs.piEN_err;
            double pEmin = fabs(obs.piEE_obs) - fISpiE * obs.piEE_err;
            if (pNmin < 0) pNmin = 0;
            if (pEmin < 0) pEmin = 0;
            double pNmax = fabs(obs.piEN_obs) + fISpiE * obs.piEN_err;
            double pEmax = fabs(obs.piEE_obs) + fISpiE * obs.piEE_err;
            obs.piE_min = sqrt(pNmin*pNmin + pEmin*pEmin);
            obs.piE_max = sqrt(pNmax*pNmax + pEmax*pEmax);
            if (obs.piE_min <= 0) obs.piE_min = 1e-10;
        }
    }

    return obs;
}

} // namespace genulens
