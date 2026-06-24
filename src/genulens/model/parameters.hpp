#pragma once

namespace genulens::model {

struct IMFParameters {
    double m0 = 1.0;
    double m1 = 0.859770466578045;
    double m2 = 0.08;
    double m3 = 0.01;
    double mbr = 0.001;
    double ml = 0.001;
    double mu = 120.0;
    double alpha0 = -2.32279457078378;
    double alpha1 = -2.32279457078378;
    double alpha2 = -1.13449983242887;
    double alpha3 = -0.175862190587576;
    double alpha4 = -0.175862190587576;
    double alpha5 = -0.175862190587576;
};

struct DensityParameters {
    int disk = 2;
    double rho_t0 = 0.042;
    int h_disk = 0;
    int add_x = 5;
    int model = 5;
    double r0 = 8160.0;
    double theta_d = 27.0;
    double frho0b = 0.839014514507754;
    double rc = 2631.78535429573;
    double zb_c = 1e+6;

    double x0 = 930.623146993329;
    double y0 = 370.784386649364;
    double z0 = 239.547516030578;
    double c1 = 1.20011972384328;
    double c2 = 4.09326795684828;
    double c3 = 1.0000;

    double x0_x = 278.027059842233;
    double y0_x = 176.318528789193;
    double z0_x = 286.791941602401;
    double c1_x = 1.3087131258784;
    double c2_x = 2.21745322869032;
    double b_zx = 1.37774815817195;
    double f_x = 1.43975636704683;
    double rc_x = 1301.63829617294;
    double b_zy = 0.0;

    int stellar_halo = 1;
    double rho0_sh_ms = 9.32e-06;
};

struct KinematicsParameters {
    double omega_p = 47.4105844018699;
    int model_vb = 5;
    double x0_vb = 858.106595717275;
    double y0_vb = 3217.04987721548;
    double z0_vb = 950.690583433628;
    double c1_vb = 4.25236641149869;
    double c2_vb = 1.02531652066343;
    double c3_vb = 1.0;
    double sigx_vb = 151.854794853683;
    double sigy_vb = 78.0278905748233;
    double sigz_vb = 81.9641955092164;
    double sigx_vb0 = 63.9939241108675;
    double sigy_vb0 = 75.8180486866697;
    double sigz_vb0 = 71.2336430487113;
    double vx_stream = 43.0364707040617;
    double y0_stream = 406.558313420815;

    int model_vbz = 5;
    double x0_vbz = 558.430182718529;
    double y0_vbz = 2003.21703656302;
    double z0_vbz = 3823.20855045157;
    double c1_vbz = 3.71001266000693;
    double c2_vbz = 1.07455173734341;
    double c3_vbz = 1.0;

    double hsig_ut = 14300.0;
    double hsig_wt = 5900.0;
    double hsig_uT = 180000.0;
    double hsig_wT = 9400.0;
    double beta_u = 0.32;
    double beta_w = 0.77;
    double sig_u10d = 42.0;
    double sig_w10d = 24.4;
    double sig_u0td = 75.0;
    double sig_w0td = 49.2;

    double sig_u_sh = 131.0;
    double sig_v_sh = 106.0;
    double sig_w_sh = 85.0;
};

struct NsdParameters {
    int enabled = -1;
    int x0 = 250;
    int y0 = 125;
    int z0 = 50;
    double mass = 0.0;
};

struct BhKickParameters {
    int mix_disk_kick = 0;
    double kick_bh = 100.0;
    double kick_ns = 350.0;
    int disk_scale_height = 0;
    int bar_scale_height = 0;
    int fix_disk_scale_length = 0;
    double disk_scale_length = 9660.0;
    double beta = 0.820;
    int use_sigma_correction = 0;
};

struct ModelParameters {
    IMFParameters imf;
    DensityParameters density;
    KinematicsParameters kinematics;
    NsdParameters nsd;
    BhKickParameters bh_kick;
};

const ModelParameters &default_model_parameters();

} // namespace genulens::model
