#include "genulens/simulation/initialize.hpp"

#include "genulens/cli/option.h"
#include "genulens/constants.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/rng.hpp"

#include <cmath>
#include <memory>

namespace genulens {

RunContext Initializer::create_context() const
{
    return {};
}

void Initializer::initialize_rng(RunContext &context, int argc, char **argv) const
{
    context.seed = getOptioni(argc, argv, "seed", 1, 12304357);
    context.runtime.rng = std::make_unique<RandomEngine>(static_cast<unsigned long>(context.seed));
}

void Initializer::read_model_options(RunContext &context, int argc, char **argv) const
{
    auto &density = context.density;
    auto &kinematics = context.kinematics;
    auto &imf = context.imf_options;
    const auto &default_imf = model::default_model_parameters().imf;

    imf.m0 = getOptiond(argc, argv, "M0", 1, default_imf.m0);
    imf.m1 = getOptiond(argc, argv, "M1", 1, default_imf.m1);
    imf.m2 = getOptiond(argc, argv, "M2", 1, default_imf.m2);
    imf.m3 = getOptiond(argc, argv, "M3", 1, default_imf.m3);
    imf.ml = getOptiond(argc, argv, "Ml", 1, default_imf.ml);
    imf.mu = getOptiond(argc, argv, "Mu", 1, default_imf.mu);
    imf.alpha1 = getOptiond(argc, argv, "alpha1", 1, default_imf.alpha1);
    imf.alpha2 = getOptiond(argc, argv, "alpha2", 1, default_imf.alpha2);
    imf.alpha3 = getOptiond(argc, argv, "alpha3", 1, default_imf.alpha3);
    imf.alpha0 = getOptiond(argc, argv, "alpha0", 1, default_imf.alpha0);
    imf.alpha4 = getOptiond(argc, argv, "alpha4", 1, default_imf.alpha4);

    density.DISK = getOptiond(argc, argv, "DISK", 1, 2);
    density.rhot0 = getOptiond(argc, argv, "rhot0", 1, 0.042);
    density.hDISK = getOptiond(argc, argv, "hDISK", 1, 0);
    density.addX = getOptiond(argc, argv, "addX", 1, 5);
    density.model = getOptiond(argc, argv, "model", 1, 5);
    density.R0 = getOptiond(argc, argv, "R0", 1, 8160);
    density.thetaD = getOptiond(argc, argv, "thetaD", 1, 27);
    density.frho0b = getOptiond(argc, argv, "frho0b", 1, 0.839014514507754);
    density.Rc = getOptiond(argc, argv, "Rc", 1, 2631.78535429573);
    density.zb_c = getOptiond(argc, argv, "zb_c", 1, 1e+6);
    if (density.model >= 4 && density.model <= 8) {
        density.x0_1 = getOptiond(argc, argv, "x0", 1, 930.623146993329);
        density.y0_1 = getOptiond(argc, argv, "y0", 1, 370.784386649364);
        density.z0_1 = getOptiond(argc, argv, "z0", 1, 239.547516030578);
        density.C1 = getOptiond(argc, argv, "C1", 1, 1.20011972384328);
        density.C2 = getOptiond(argc, argv, "C2", 1, 4.09326795684828);
        density.C3 = getOptiond(argc, argv, "C3", 1, 1.0000);
    }
    if (density.addX >= 5) {
        density.x0_X = getOptiond(argc, argv, "x0_X", 1, 278.027059842233);
        density.y0_X = getOptiond(argc, argv, "y0_X", 1, 176.318528789193);
        density.z0_X = getOptiond(argc, argv, "z0_X", 1, 286.791941602401);
        density.C1_X = getOptiond(argc, argv, "C1_X", 1, 1.3087131258784);
        density.C2_X = getOptiond(argc, argv, "C2_X", 1, 2.21745322869032);
        density.b_zX = getOptiond(argc, argv, "b_zX", 1, 1.37774815817195);
        density.fX = getOptiond(argc, argv, "fX", 1, 1.43975636704683);
        density.Rc_X = getOptiond(argc, argv, "Rc_X", 1, 1301.63829617294);
    }
    density.b_zY = getOptiond(argc, argv, "b_zY", 1, 0);

    kinematics.Omega_p = getOptiond(argc, argv, "Omega_p", 1, 47.4105844018699);
    kinematics.model_vb = getOptiond(argc, argv, "model_vb", 1, 5);
    kinematics.x0_vb = getOptiond(argc, argv, "x0_vb", 1, 858.106595717275);
    kinematics.y0_vb = getOptiond(argc, argv, "y0_vb", 1, 3217.04987721548);
    kinematics.z0_vb = getOptiond(argc, argv, "z0_vb", 1, 950.690583433628);
    kinematics.C1_vb = getOptiond(argc, argv, "C1_vb", 1, 4.25236641149869);
    kinematics.C2_vb = getOptiond(argc, argv, "C2_vb", 1, 1.02531652066343);
    kinematics.C3_vb = getOptiond(argc, argv, "C3_vb", 1, 1);
    kinematics.sigx_vb = getOptiond(argc, argv, "sigx_vb", 1, 151.854794853683);
    kinematics.sigy_vb = getOptiond(argc, argv, "sigy_vb", 1, 78.0278905748233);
    kinematics.sigz_vb = getOptiond(argc, argv, "sigz_vb", 1, 81.9641955092164);
    kinematics.sigx_vb0 = getOptiond(argc, argv, "sigx_vb0", 1, 63.9939241108675);
    kinematics.sigy_vb0 = getOptiond(argc, argv, "sigy_vb0", 1, 75.8180486866697);
    kinematics.sigz_vb0 = getOptiond(argc, argv, "sigz_vb0", 1, 71.2336430487113);
    kinematics.vx_str = getOptiond(argc, argv, "vx_str", 1, 43.0364707040617);
    kinematics.y0_str = getOptiond(argc, argv, "y0_str", 1, 406.558313420815);
    kinematics.model_vbz = getOptiond(argc, argv, "model_vbz", 1, 5);
    kinematics.x0_vbz = getOptiond(argc, argv, "x0_vbz", 1, 558.430182718529);
    kinematics.y0_vbz = getOptiond(argc, argv, "y0_vbz", 1, 2003.21703656302);
    kinematics.z0_vbz = getOptiond(argc, argv, "z0_vbz", 1, 3823.20855045157);
    kinematics.C1_vbz = getOptiond(argc, argv, "C1_vbz", 1, 3.71001266000693);
    kinematics.C2_vbz = getOptiond(argc, argv, "C2_vbz", 1, 1.07455173734341);
    kinematics.C3_vbz = getOptiond(argc, argv, "C3_vbz", 1, 1);

    kinematics.hsigUt = getOptiond(argc, argv, "hsigUt", 1, 14300);
    kinematics.hsigWt = getOptiond(argc, argv, "hsigWt", 1, 5900);
    kinematics.hsigUT = getOptiond(argc, argv, "hsigUT", 1, 180000);
    kinematics.hsigWT = getOptiond(argc, argv, "hsigWT", 1, 9400);
    kinematics.betaU = getOptiond(argc, argv, "betaU", 1, 0.32);
    kinematics.betaW = getOptiond(argc, argv, "betaW", 1, 0.77);
    kinematics.sigU10d = getOptiond(argc, argv, "sigU10d", 1, 42.0);
    kinematics.sigW10d = getOptiond(argc, argv, "sigW10d", 1, 24.4);
    kinematics.sigU0td = getOptiond(argc, argv, "sigU0td", 1, 75.0);
    kinematics.sigW0td = getOptiond(argc, argv, "sigW0td", 1, 49.2);

    context.B14disk = getOptiond(argc, argv, "B14disk", 1, 0);
    const int B14bar = getOptiond(argc, argv, "B14bar", 1, 0);
    const int E_fg0 = getOptiond(argc, argv, "E_fg0", 1, 0);
    const int G_fg0 = getOptiond(argc, argv, "G_fg0", 1, 0);
    const int EXE_fg0 = getOptiond(argc, argv, "EXE_fg0", 1, 0);
    const int GXG_fg0 = getOptiond(argc, argv, "GXG_fg0", 1, 0);
    if (context.B14disk == 1) {
        kinematics.vxsun = -12.7, kinematics.vysun = 218.0 + 24.0, kinematics.vzsun = 7.25;
        density.DISK = 1;
    }
    if (B14bar == 1) {
        imf.m0 = 1.0, imf.m1 = 0.7, imf.m2 = 0.08, imf.m3 = 0.01;
        imf.alpha1 = -2.0, imf.alpha2 = -1.3, imf.alpha3 = -0.5;
        imf.alpha0 = imf.alpha1, imf.alpha4 = imf.alpha3;
        density.model = 6, density.addX = 0, context.B14vbar = 1;
        density.R0 = 8200, density.thetaD = 20, density.x0_1 = 1580.0, density.y0_1 = 620.0;
        density.z0_1 = 430.0, density.Rc = 2400.0, density.C1 = 2, density.C2 = 4, density.C3 = 1;
        density.frho0b = 1.173;
        kinematics.Omega_p = 50, kinematics.sigx_vb = 114.0, kinematics.sigy_vb = 103.8,
        kinematics.sigz_vb = 96.4;
        kinematics.x0_vb = kinematics.y0_vb = kinematics.z0_vb = 500000;
        kinematics.x0_vbz = kinematics.y0_vbz = kinematics.z0_vbz = 500000;
    }
    if (E_fg0 == 1) {
        density.model = 5, density.addX = 0;
        imf.m0 = 1.0, imf.m1 = 0.843651488650385, imf.m2 = 0.08, imf.m3 = 0.01;
        imf.alpha1 = -2.30708461042964, imf.alpha2 = -1.09811414023325, imf.alpha3 = -0.176687444667866;
        imf.alpha0 = imf.alpha1, imf.alpha4 = imf.alpha3;
        density.R0 = 8160, density.thetaD = 27, density.frho0b = 0.847695765083198,
        density.Rc = 2804.94024639663;
        density.x0_1 = 668.323640191308, density.y0_1 = 277.674592258175,
        density.z0_1 = 235.344943180979, density.C1 = 1.40903573470129,
        density.C2 = 3.3497118832179, density.C3 = 1;
        kinematics.model_vb = 5, kinematics.model_vbz = 5;
        kinematics.Omega_p = 49.5149910609312, kinematics.vx_str = 48.7482280102778,
        kinematics.y0_str = 392.515724264323, kinematics.sigx_vb = 156.055410564041,
        kinematics.sigy_vb = 83.8197043324931, kinematics.sigz_vb = 86.3564038759999,
        kinematics.sigx_vb0 = 63.8292191277825, kinematics.sigy_vb0 = 74.9469462226124,
        kinematics.sigz_vb0 = 72.3085487545662, kinematics.x0_vb = 823.387929122523,
        kinematics.y0_vb = 9288.51482678556, kinematics.z0_vb = 864.479916419292,
        kinematics.C1_vb = 3.82820123451928, kinematics.C2_vb = 1.00573720627546,
        kinematics.x0_vbz = 511.063328964278, kinematics.y0_vbz = 2896.01606378595,
        kinematics.z0_vbz = 2189.7664883434, kinematics.C1_vbz = 3.04214421342047,
        kinematics.C2_vbz = 1.00609904766722;
    }
    if (G_fg0 == 1) {
        density.model = 6, density.addX = 0;
        imf.m0 = 1.0, imf.m1 = 0.896557393600988, imf.m2 = 0.08, imf.m3 = 0.01;
        imf.alpha1 = -2.39628188518525, imf.alpha2 = -1.18451896148506, imf.alpha3 = 0.168672130848533;
        imf.alpha0 = imf.alpha1, imf.alpha4 = imf.alpha3;
        density.R0 = 8160, density.thetaD = 27, density.frho0b = 0.777347874844233,
        density.Rc = 4838.85613149588;
        density.x0_1 = 1025.42128394916, density.y0_1 = 457.419718281149,
        density.z0_1 = 396.048253079423, density.C1 = 2.00928445577057,
        density.C2 = 3.9678518191928, density.C3 = 1;
        kinematics.model_vb = 5, kinematics.model_vbz = 5;
        kinematics.Omega_p = 40.5174879673548, kinematics.vx_str = 11.9026090372449,
        kinematics.y0_str = 20.1384817812277, kinematics.sigx_vb = 136.435675357212,
        kinematics.sigy_vb = 109.313291840218, kinematics.sigz_vb = 101.291432907346,
        kinematics.sigx_vb0 = 76.0453005937702, kinematics.sigy_vb0 = 67.9783132842431,
        kinematics.sigz_vb0 = 74.7117386554542, kinematics.x0_vb = 1031.18302251324,
        kinematics.y0_vb = 2145.45565210108, kinematics.z0_vb = 727.233943973984,
        kinematics.C1_vb = 4.9302429910108, kinematics.C2_vb = 1.04038121792228,
        kinematics.x0_vbz = 517.854475368706, kinematics.y0_vbz = 1436.21008855387,
        kinematics.z0_vbz = 1095.79181359292, kinematics.C1_vbz = 2.3091601785779,
        kinematics.C2_vbz = 1.03832670354301;
    }
    if (EXE_fg0 == 1) {
        density.model = 5, density.addX = 5;
        imf.m0 = 1.0, imf.m1 = 0.859770466578045, imf.m2 = 0.08, imf.m3 = 0.01;
        imf.alpha1 = -2.32279457078378, imf.alpha2 = -1.13449983242887, imf.alpha3 = -0.175862190587576;
        imf.alpha0 = imf.alpha1, imf.alpha4 = imf.alpha3;
        density.R0 = 8160, density.thetaD = 27, density.frho0b = 0.839014514507754,
        density.Rc = 2631.78535429573;
        density.x0_1 = 930.623146993329, density.y0_1 = 370.784386649364,
        density.z0_1 = 239.547516030578, density.C1 = 1.20011972384328,
        density.C2 = 4.09326795684828, density.C3 = 1;
        kinematics.model_vb = 5, kinematics.model_vbz = 5;
        kinematics.Omega_p = 47.4105844018699, kinematics.vx_str = 43.0364707040617,
        kinematics.y0_str = 406.558313420815, kinematics.sigx_vb = 151.854794853683,
        kinematics.sigy_vb = 78.0278905748233, kinematics.sigz_vb = 81.9641955092164,
        kinematics.sigx_vb0 = 63.9939241108675, kinematics.sigy_vb0 = 75.8180486866697,
        kinematics.sigz_vb0 = 71.2336430487113, kinematics.x0_vb = 858.106595717275,
        kinematics.y0_vb = 3217.04987721548, kinematics.z0_vb = 950.690583433628,
        kinematics.C1_vb = 4.25236641149869, kinematics.C2_vb = 1.02531652066343,
        kinematics.x0_vbz = 558.430182718529, kinematics.y0_vbz = 2003.21703656302,
        kinematics.z0_vbz = 3823.20855045157, kinematics.C1_vbz = 3.71001266000693,
        kinematics.C2_vbz = 1.07455173734341;
        density.x0_X = 278.027059842233, density.y0_X = 176.318528789193,
        density.z0_X = 286.791941602401, density.C1_X = 1.3087131258784,
        density.C2_X = 2.21745322869032, density.b_zX = 1.37774815817195,
        density.fX = 1.43975636704683, density.Rc_X = 1301.63829617294;
    }
    if (GXG_fg0 == 1) {
        density.model = 6, density.addX = 6;
        imf.m0 = 1.0, imf.m1 = 0.901747918318042, imf.m2 = 0.08, imf.m3 = 0.01;
        imf.alpha1 = -2.32055781291126, imf.alpha2 = -1.16146692073597, imf.alpha3 = -0.222751835826612;
        imf.alpha0 = imf.alpha1, imf.alpha4 = imf.alpha3;
        density.R0 = 8160, density.thetaD = 27, density.frho0b = 0.861982105059042,
        density.Rc = 2834.43172768484;
        density.x0_1 = 1564.78976595399, density.y0_1 = 721.729645984158,
        density.z0_1 = 494.669973292979, density.C1 = 1.20141097225,
        density.C2 = 3.09254667088709, density.C3 = 1;
        kinematics.model_vb = 5, kinematics.model_vbz = 5;
        kinematics.Omega_p = 45.9061365175252, kinematics.vx_str = 28.250608437116,
        kinematics.y0_str = 11.4387290790323, kinematics.sigx_vb = 154.984185643613,
        kinematics.sigy_vb = 78.4783157632334, kinematics.sigz_vb = 83.2424209150283,
        kinematics.sigx_vb0 = 63.3834790223473, kinematics.sigy_vb0 = 75.1951371572303,
        kinematics.sigz_vb0 = 69.6076680158332, kinematics.x0_vb = 939.470002303028,
        kinematics.y0_vb = 4228.61947632437, kinematics.z0_vb = 883.716365308057,
        kinematics.C1_vb = 4.59067123072475, kinematics.C2_vb = 1.00961963171066,
        kinematics.x0_vbz = 699.073733500672, kinematics.y0_vbz = 1729.91970395558,
        kinematics.z0_vbz = 2028.24030134845, kinematics.C1_vbz = 4.84589813971794,
        kinematics.C2_vbz = 1.01718557457505;
        density.x0_X = 755.975821023038, density.y0_X = 312.17136920671,
        density.z0_X = 399.287597819655, density.C1_X = 1.21131134854495,
        density.C2_X = 1.30388556329566, density.b_zX = 1.37711800325276,
        density.fX = 2.99985800759016, density.Rc_X = 5174.00544959931;
    }
}

void Initializer::finalize_spatial_model(RunContext &context, int argc, char **argv) const
{
    auto &density = context.density;
    auto &kinematics = context.kinematics;

    density.costheta = std::cos(density.thetaD / 180.0 * kPi);
    density.sintheta = std::sin(density.thetaD / 180.0 * kPi);

    auto &spatial = context.spatial;
    spatial.center_on_sgr_a = getOptioni(argc, argv, "CenSgrA", 1, spatial.center_on_sgr_a);
    if (spatial.center_on_sgr_a == 1) {
        context.xyzSgrA = model::CoordinateTransformer().distance_l_b_to_xyz(
            density.R0, spatial.l_sgr_a, spatial.b_sgr_a, density.R0);
    }

    density.SH = getOptiond(argc, argv, "SH", 1, density.SH);
    density.rho0SHMS = getOptiond(argc, argv, "rho0SHMS", 1, 9.32e-06);
    kinematics.sigU_SH = getOptiond(argc, argv, "sigU_SH", 1, 131);
    kinematics.sigV_SH = getOptiond(argc, argv, "sigV_SH", 1, 106);
    kinematics.sigW_SH = getOptiond(argc, argv, "sigW_SH", 1, 85);
    if (density.SH == 0) density.rho0SHMS = 0.0;
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
