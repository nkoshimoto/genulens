#pragma once

namespace genulens {

struct ObservationConfig {
    double tE_obs = 54.5, tE_err = 99999999999.0, fe_tE = 0;
    int tE_det = 0;
    double tE_min = 0, tE_max = 0;

    double thetaE_obs = 0, thetaE_err = 0, fe_thetaE = 0;
    int thetaE_det = 0;
    double thetaE_min = 0, thetaE_max = 0;

    double piE_obs = 0, piE_err = 0, fe_piE = 0;
    int piE_det = 0;
    double piE_min = 0, piE_max = 0;

    double piEN_obs = 0, piEN_err = 0, fe_piEN = 0;
    double piEE_obs = 0, piEE_err = 0, fe_piEE = 0;

    double musl_obs = 0, musl_err = 0, fe_musl = 0;
    double musb_obs = 0, musb_err = 0, fe_musb = 0;
    double musN_obs = 0, musN_err = 0, fe_musN = 0;
    double musE_obs = 0, musE_err = 0, fe_musE = 0;
    int musRCG = 0;

    double muhelN_obs = 0, muhelN_err = 0, fe_muhelN = 0;
    double muhelE_obs = 0, muhelE_err = 0, fe_muhelE = 0;

    double IL_obs = 14.0, IL_err = 0.01, fe_IL = 0;
    int IL_det = 2;
    double KL_obs = 0, KL_err = 0, fe_KL = 0;
    int KL_det = 0;

    double u0_obs = 0;

    static ObservationConfig from_cli(int argc, char **argv, int uniform_likelihood);
};

} // namespace genulens
