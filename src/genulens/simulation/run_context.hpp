#pragma once

#include "genulens/rng.hpp"

#include <array>
#include <memory>

namespace genulens {

struct ScientificRuntimeState {
    std::unique_ptr<RandomEngine> rng;
};

struct SamplingOptions {
    long n_simu = 100000;
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

struct StellarPopulationState {
    double tSFR = 7.0;
    double MiniWDmax = 9.0;

    std::array<int, 250> agesD{};
    std::array<int, 50> agesB{};
    std::array<int, 10> agesND{};
    std::array<double, 250> MinidieD{};
    std::array<double, 50> MinidieB{};
    std::array<double, 10> MinidieND{};
    int nageD = 0;
    int nageB = 0;
    int nageND = 0;
    double mageB = 9.0;
    double sageB = 1.0;
    double mageND = 7.0;
    double sageND = 1.0;

    int nm = 0;
    double logMst = 0.0;
    double dlogM = 0.0;
};

struct DensityState {
    int ncomp = 11;
    double rhot0 = 0.0;

    double fb_MS = 1.62 / 2.07;
    double m2nb_MS = 1.0 / 0.227943;
    double m2nb_WD = 1.0 / 0.847318;
    double nMS2nRGb = 2.33232e-03;
    double rho0b = 0.0;
    double n0MSb = 0.0;
    double n0RGb = 0.0;
    double n0b = 0.0;

    int ND = 0;
    int x0ND = 250;
    int y0ND = 125;
    int z0ND = 50;
    double C1ND = 2.0;
    double rho0ND = 0.0;
    double n0MSND = 0.0;
    double n0RGND = 0.0;
    double n0ND = 0.0;
    double fND_MS = 0.0;
    double m2nND_MS = 0.0;
    double m2nND_WD = 0.0;
    double nMS2nRGND = 0.0;

    int SH = 1;
    double rho0SHMS = 9.32e-06;
    double epsSH = 0.76;
    double alphaSH = 2.44;
    double acSH2 = 500.0 * 500.0;
    double rho0SH = 0.0;
    double n0MSSH = 0.0;
    double n0RGSH = 0.0;
    double n0SH = 0.0;

    std::array<double, 8> rho0d{{5.16e-03 + 3.10e-04, 5.00e-03 + 5.09e-04, 3.85e-03 + 5.42e-04, 3.18e-03 + 5.54e-04,
                                 5.84e-03 + 1.21e-03, 6.24e-03 + 1.51e-03, 1.27e-02 + 3.49e-03, 1.68e-03 + 6.02e-04}};
    std::array<double, 8> n0d{{1.51e-02 + 1.12e-04, 1.66e-02 + 3.22e-04, 1.40e-02 + 4.39e-04, 1.22e-02 + 5.15e-04,
                               2.36e-02 + 1.25e-03, 2.63e-02 + 1.67e-03, 5.55e-02 + 4.08e-03, 7.91e-03 + 7.81e-04}};
    std::array<double, 8> n0MSd{{1.51e-02, 1.66e-02, 1.40e-02, 1.22e-02, 2.36e-02, 2.63e-02, 5.55e-02, 7.91e-03}};
    std::array<double, 8> n0RGd{{7.09e-06, 3.40e-05, 4.32e-05, 2.16e-05, 6.60e-05, 6.19e-05, 1.29e-04, 9.38e-06}};

    std::array<double, 3> y0d{};
    std::array<int, 3> Rd{{5000, 2600, 2200}};
    int Rh = 3740;
    int Rdbreak = 5300;
    int nh = 1;
    std::array<double, 8> zd{{61.47, 141.84, 224.26, 292.36, 372.85, 440.71, 445.37, 903.12}};
    std::array<double, 8> zd45{{36.88, 85.10, 134.55, 175.41, 223.71, 264.42, 267.22, 903.12}};

    int DISK = 0;
    int hDISK = 0;
    int addX = 0;
    int model = 0;
    double R0 = 0.0;
    double thetaD = 0.0;
    double x0_1 = 0.0;
    double y0_1 = 0.0;
    double z0_1 = 0.0;
    double C1 = 0.0;
    double C2 = 0.0;
    double C3 = 0.0;
    double Rc = 0.0;
    double frho0b = 0.0;
    double costheta = 0.0;
    double sintheta = 0.0;
    double zb_c = 0.0;
    double x0_X = 0.0;
    double y0_X = 0.0;
    double z0_X = 0.0;
    double C1_X = 0.0;
    double C2_X = 0.0;
    double b_zX = 0.0;
    double fX = 0.0;
    double Rsin = 0.0;
    double b_zY = 0.0;
    double Rc_X = 0.0;

    double *lDs = nullptr;
    double *bDs = nullptr;
};

struct LuminosityFunctionState {
    int nMIs = 0;
    int nVIs = 0;
    double *MIs = nullptr;
    double **CumuN_MIs = nullptr;
    double dILF = 0.0;
    double *VIs = nullptr;
    double ***f_VI_Is = nullptr;
    double dVILF = 0.0;
};

struct KinematicsState {
    int nVcs = 0;
    std::array<double, 60> Rcs{};
    std::array<double, 60> Vcs{};

    double vxsun = -10.0;
    double Vsun = 11.0;
    double vzsun = 7.0;
    double vysun = 243.0;

    double ****fgsShu = nullptr;
    double ****PRRgShus = nullptr;
    double ****cumu_PRRgs = nullptr;
    int ***n_fgsShu = nullptr;
    int ****kptiles = nullptr;
    double hsigUt = 0.0;
    double hsigWt = 0.0;
    double hsigUT = 0.0;
    double hsigWT = 0.0;
    double betaU = 0.0;
    double betaW = 0.0;
    double sigU10d = 0.0;
    double sigW10d = 0.0;
    double sigU0td = 0.0;
    double sigW0td = 0.0;
    std::array<double, 8> medtauds{{0.075273, 0.586449, 1.516357, 2.516884, 4.068387, 6.069263, 8.656024, 12.0}};

    int zstShu = 0;
    int zenShu = 3600;
    int dzShu = 200;
    int RstShu = 500;
    int RenShu = 12200;
    int dRShu = 100;

    int model_vb = 0;
    int model_vbz = 0;
    double Omega_p = 0.0;
    double x0_vb = 0.0;
    double y0_vb = 0.0;
    double z0_vb = 0.0;
    double C1_vb = 0.0;
    double C2_vb = 0.0;
    double C3_vb = 0.0;
    double sigx_vb = 0.0;
    double sigy_vb = 0.0;
    double sigz_vb = 0.0;
    double vx_str = 0.0;
    double y0_str = 0.0;
    double sigx_vb0 = 0.0;
    double sigy_vb0 = 0.0;
    double sigz_vb0 = 0.0;
    double x0_vbz = 0.0;
    double y0_vbz = 0.0;
    double z0_vbz = 0.0;
    double C1_vbz = 0.0;
    double C2_vbz = 0.0;
    double C3_vbz = 0.0;

    double sigU_SH = 0.0;
    double sigV_SH = 0.0;
    double sigW_SH = 0.0;
};

struct NSDMomentsState {
    double **logrhoNDs = nullptr;
    double **vphiNDs = nullptr;
    double ***logsigvNDs = nullptr;
    double **corRzNDs = nullptr;
    double zstND = 0.0;
    double zenND = 400.0;
    double dzND = 5.0;
    double RstND = 0.0;
    double RenND = 1000.0;
    double dRND = 5.0;
    int nzND = 0;
    int nRND = 0;
};

struct GenulensRunContext {
    ScientificRuntimeState runtime;
    StellarPopulationState stellar;
    DensityState density;
    LuminosityFunctionState luminosity;
    KinematicsState kinematics;
    NSDMomentsState nsd_moments;

    long seed = 0;
    int B14disk = 0;
    int B14vbar = 0;
    std::array<double, 3> xyzSgrA{};
    SamplingOptions sampling;
};

} // namespace genulens
