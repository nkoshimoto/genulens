/* Generate microlensing events following the Galactic model developed by Koshimoto, Baba & Bennett (2021).
 * N. Koshimoto wrote the original .c version and C. Ranc converted it into .cpp.
 * This is version 1.2.1+ of genulens (refactored). */
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdexcept>
#include <utility>
#include "genulens/cli/option.h"
#include "genulens/io/input_data.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/extinction.hpp"
#include "genulens/model/forward_source.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/rng.hpp"
#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/los_density_grid.hpp"
#include "genulens/simulation/mass_function.hpp"
#include "genulens/simulation/observation_config.hpp"
#include "genulens/simulation/prepared_simulation.hpp"
#include "genulens/simulation/sampler.hpp"
#include "genulens/simulation/velocity_distribution.hpp"

#define fopen(path, mode) genulens::open_input_file((path), (mode))

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"

#define PI 3.1415926535897932385

namespace genulens {

namespace {

gmodel::ForwardSourceGenerator make_forward_source_generator(const GenulensConfig &config)
{
    if (config.forward_source.photometry == "prime") {
        return gmodel::ForwardSourceGenerator::load_default_prime(config.model.imf);
    }
    if (config.forward_source.photometry == "roman") {
        return gmodel::ForwardSourceGenerator::load_default_roman(config.model.imf);
    }
    throw std::runtime_error("unknown forward source photometry: " + config.forward_source.photometry);
}

std::vector<std::string> forward_source_bands_for_config(const GenulensConfig &config)
{
    if (config.forward_source.enabled == 0) return {};
    return make_forward_source_generator(config).bands();
}

void apply_typed_model_parameters(RunContext &context, const model::ModelParameters &model)
{
    const auto &density_cfg = model.density;
    auto &density = context.density;
    density.DISK = density_cfg.disk;
    density.rhot0 = density_cfg.rho_t0;
    density.hDISK = density_cfg.h_disk;
    density.addX = density_cfg.add_x;
    density.model = density_cfg.model;
    density.R0 = density_cfg.r0;
    density.thetaD = density_cfg.theta_d;
    density.frho0b = density_cfg.frho0b;
    density.Rc = density_cfg.rc;
    density.zb_c = density_cfg.zb_c;
    density.x0_1 = density_cfg.x0;
    density.y0_1 = density_cfg.y0;
    density.z0_1 = density_cfg.z0;
    density.C1 = density_cfg.c1;
    density.C2 = density_cfg.c2;
    density.C3 = density_cfg.c3;
    density.x0_X = density_cfg.x0_x;
    density.y0_X = density_cfg.y0_x;
    density.z0_X = density_cfg.z0_x;
    density.C1_X = density_cfg.c1_x;
    density.C2_X = density_cfg.c2_x;
    density.b_zX = density_cfg.b_zx;
    density.fX = density_cfg.f_x;
    density.Rc_X = density_cfg.rc_x;
    density.b_zY = density_cfg.b_zy;
    density.SH = density_cfg.stellar_halo;
    density.rho0SHMS = density_cfg.rho0_sh_ms;
    if (density.SH == 0) density.rho0SHMS = 0.0;

    const auto &kinematics_cfg = model.kinematics;
    auto &kinematics = context.kinematics;
    kinematics.Omega_p = kinematics_cfg.omega_p;
    kinematics.model_vb = kinematics_cfg.model_vb;
    kinematics.x0_vb = kinematics_cfg.x0_vb;
    kinematics.y0_vb = kinematics_cfg.y0_vb;
    kinematics.z0_vb = kinematics_cfg.z0_vb;
    kinematics.C1_vb = kinematics_cfg.c1_vb;
    kinematics.C2_vb = kinematics_cfg.c2_vb;
    kinematics.C3_vb = kinematics_cfg.c3_vb;
    kinematics.sigx_vb = kinematics_cfg.sigx_vb;
    kinematics.sigy_vb = kinematics_cfg.sigy_vb;
    kinematics.sigz_vb = kinematics_cfg.sigz_vb;
    kinematics.sigx_vb0 = kinematics_cfg.sigx_vb0;
    kinematics.sigy_vb0 = kinematics_cfg.sigy_vb0;
    kinematics.sigz_vb0 = kinematics_cfg.sigz_vb0;
    kinematics.vx_str = kinematics_cfg.vx_stream;
    kinematics.y0_str = kinematics_cfg.y0_stream;
    kinematics.model_vbz = kinematics_cfg.model_vbz;
    kinematics.x0_vbz = kinematics_cfg.x0_vbz;
    kinematics.y0_vbz = kinematics_cfg.y0_vbz;
    kinematics.z0_vbz = kinematics_cfg.z0_vbz;
    kinematics.C1_vbz = kinematics_cfg.c1_vbz;
    kinematics.C2_vbz = kinematics_cfg.c2_vbz;
    kinematics.C3_vbz = kinematics_cfg.c3_vbz;
    kinematics.hsigUt = kinematics_cfg.hsig_ut;
    kinematics.hsigWt = kinematics_cfg.hsig_wt;
    kinematics.hsigUT = kinematics_cfg.hsig_uT;
    kinematics.hsigWT = kinematics_cfg.hsig_wT;
    kinematics.betaU = kinematics_cfg.beta_u;
    kinematics.betaW = kinematics_cfg.beta_w;
    kinematics.sigU10d = kinematics_cfg.sig_u10d;
    kinematics.sigW10d = kinematics_cfg.sig_w10d;
    kinematics.sigU0td = kinematics_cfg.sig_u0td;
    kinematics.sigW0td = kinematics_cfg.sig_w0td;
    kinematics.sigU_SH = kinematics_cfg.sig_u_sh;
    kinematics.sigV_SH = kinematics_cfg.sig_v_sh;
    kinematics.sigW_SH = kinematics_cfg.sig_w_sh;
}

} // namespace



int run_sampler_impl(RunContext &context,
                     const GenulensConfig *typed_config,
                     int argc, char **argv,
                     LikelihoodFunction custom_likelihood,
                     EventSampler::EventSink event_sink,
                     bool emit_cli_output)
{
    char curdir[200];
    getcwd(curdir, 200);
    PreparedSimulation prepared;

    // ------- Model options -------
    int CheckD = getOptiond(argc, argv, "CheckD", 1, 0);
    Initializer().read_model_options(context, argc, argv);
    if (typed_config) {
        const auto &imf = typed_config->model.imf;
        context.imf_options.m0 = imf.m0;
        context.imf_options.m1 = imf.m1;
        context.imf_options.m2 = imf.m2;
        context.imf_options.m3 = imf.m3;
        context.imf_options.ml = imf.ml;
        context.imf_options.mu = imf.mu;
        context.imf_options.alpha0 = imf.alpha0;
        context.imf_options.alpha1 = imf.alpha1;
        context.imf_options.alpha2 = imf.alpha2;
        context.imf_options.alpha3 = imf.alpha3;
        context.imf_options.alpha4 = imf.alpha4;
        apply_typed_model_parameters(context, typed_config->model);
    }
    double M0_B = context.imf_options.m0, M1_B = context.imf_options.m1;
    double M2_B = context.imf_options.m2, M3_B = context.imf_options.m3;
    double Ml   = context.imf_options.ml, Mu   = context.imf_options.mu;
    double alpha0_B = context.imf_options.alpha0, alpha1_B = context.imf_options.alpha1;
    double alpha2_B = context.imf_options.alpha2, alpha3_B = context.imf_options.alpha3;
    double alpha4_B = context.imf_options.alpha4;
    Initializer().finalize_spatial_model(context, argc, argv);
    if (typed_config) {
        apply_typed_model_parameters(context, typed_config->model);
        context.density.costheta = std::cos(context.density.thetaD / 180.0 * PI);
        context.density.sintheta = std::sin(context.density.thetaD / 180.0 * PI);
    }

    // ------- Sightline & source-selection -------
    double lSIMU = getOptiond(argc, argv, "l",        1,  1.0);
    double bSIMU = getOptiond(argc, argv, "b",        1, -3.9);
    if (typed_config) {
        lSIMU = typed_config->l;
        bSIMU = typed_config->b;
    }
    double Isst  = getOptiond(argc, argv, "Isrange",  1, 14.0);
    double Isen  = getOptiond(argc, argv, "Isrange",  2, 21.0);
    double VIsst = getOptiond(argc, argv, "VIsrange", 1,  0.0);
    double VIsen = getOptiond(argc, argv, "VIsrange", 2,  0.0);
    double AIrc  = getOptiond(argc, argv, "AIrc",     1,  0);
    double EVIrc = getOptiond(argc, argv, "EVIrc",    1,  0);
    double DMrc  = getOptiond(argc, argv, "DMrc",     1,  0);
    double AKrc  = getOptiond(argc, argv, "AKrc",     1,  0);
    if (typed_config) {
        Isst = typed_config->source.i_min;
        Isen = typed_config->source.i_max;
        VIsst = typed_config->source.vi_min;
        VIsen = typed_config->source.vi_max;
        AIrc = typed_config->source.ai_rc;
        EVIrc = typed_config->source.evi_rc;
        DMrc = typed_config->source.dm_rc;
        AKrc = typed_config->source.ak_rc;
    }
    if (emit_cli_output && (fabs(lSIMU) > 10 || fabs(bSIMU) > 7 || fabs(bSIMU) < 1.5))
        printf("# WARNING: genulens is designed for |l| < ~10 and ~1.5 < |b| < ~7"
               " and (l,b)= (%.3f, %.3f) is outside the range.\n", lSIMU, bSIMU);

    // ------- Population runtime -------
    prepared.population.initialize_mass_function(context, context.imf_options);
    prepared.population_active = true;
    prepared.population.initialize_luminosity_functions(context, Isst, Isen, VIsst, VIsen, AIrc, EVIrc);
    prepared.luminosity_functions_active = true;
    prepared.population.read_empirical_mass_luminosity();
    double *logMass_B        = prepared.population.log_mass;
    double *PlogM_B          = prepared.population.mass_probability;
    double *PlogM_cum_norm_B = prepared.population.mass_cumulative;
    int    *imptiles_B       = prepared.population.mass_percentiles;

    // ------- Kinematic tables -------
    prepared.kinematic_tables.initialize_shu_distribution(context);
    prepared.kinematic_tables_active = true;

    VelocityDistribution vel_dist;
    vel_dist.bind(context);

    // ------- Disk normalization -------
    context.density.y0d.data()[0] = (context.density.DISK == 1) ? exp(-context.density.R0/context.density.Rd.data()[0] - pow((double)context.density.Rh/context.density.R0, context.density.nh)) : exp(-context.density.R0/context.density.Rd.data()[0]);
    context.density.y0d.data()[1] = (context.density.DISK == 1) ? exp(-context.density.R0/context.density.Rd.data()[1] - pow((double)context.density.Rh/context.density.R0, context.density.nh)) : exp(-context.density.R0/context.density.Rd.data()[1]);
    context.density.y0d.data()[2] = (context.density.DISK == 1) ? exp(-context.density.R0/context.density.Rd.data()[2] - pow((double)context.density.Rh/context.density.R0, context.density.nh)) : exp(-context.density.R0/context.density.Rd.data()[2]);

    // ------- BH kick options -------
    int    MXDkick  = getOptiond(argc, argv, "MXDkick",  1,     0);
    double vkickBH  = getOptiond(argc, argv, "vkickBH",  1,  100.0);
    double vkickNS  = getOptiond(argc, argv, "vkickNS",  1,  350.0);
    int    BHhd     = getOptiond(argc, argv, "BHhd",     1,     0);
    int    BHhb     = getOptiond(argc, argv, "BHhb",     1,     0);
    int    fixRhdBH = getOptioni(argc, argv, "fixRhdBH", 1,     0);
    double RhdBH0   = getOptiond(argc, argv, "RhdBH0",   1,  9660.0);
    double betaBH   = getOptiond(argc, argv, "betaBH",   1,  0.820);
    int    UseSigBH = getOptiond(argc, argv, "UseSigBH", 1,     0);
    int  printBHfac = getOptiond(argc, argv, "printBHfac", 1,   0);
    if (typed_config) {
        const auto &bh = typed_config->model.bh_kick;
        MXDkick = bh.mix_disk_kick;
        vkickBH = bh.kick_bh;
        vkickNS = bh.kick_ns;
        BHhd = bh.disk_scale_height;
        BHhb = bh.bar_scale_height;
        fixRhdBH = bh.fix_disk_scale_length;
        RhdBH0 = bh.disk_scale_length;
        betaBH = bh.beta;
        UseSigBH = bh.use_sigma_correction;
    }

    // ------- Header printing -------
    if (emit_cli_output) {
    printf("#   Output of \"./genulens ");
    for (int i = 1; i < argc; i++) { printf("%s", argv[i]); if (i < argc-1) printf(" "); }
    printf("\"\n");
    printf("#   You need to weight each line by the 1st value of each line from the left \n");
    printf("#---------- Parameters for IMF and Sun ----------\n");
    printf("#      IMF:  alpha0= %5.2f ( %.2f <M< %.2f ),\n", alpha0_B, M0_B, Mu);
    printf("#            alpha1= %5.2f ( %.2f <M< %.2f ),\n", alpha1_B, M1_B, M0_B);
    printf("#            alpha2= %5.2f ( %.2f <M< %.2f ),\n", alpha2_B, M2_B, M1_B);
    printf("#            alpha3= %5.2f ( %.2f <M< %.2f ),\n", alpha3_B, M3_B, M2_B);
    printf("#            alpha4= %5.2f ( %.5f <M< %.5f )\n",  alpha4_B, Ml, M3_B);
    printf("#         (R, z)sun= (%5.0f, %5.0f) pc\n", context.density.R0, 25.0);
    printf("#   (vx, vy, vz)sun= (%5.1f, %5.1f, %4.1f) km/s\n", context.kinematics.vxsun, context.kinematics.vysun, context.kinematics.vzsun);
    printf("#------------ Disk model: (DISK, hDISK, tSFR)= ( %d , %d , %.1f Gyr ) -----\n",
           context.density.DISK, context.density.hDISK, context.stellar.tSFR);
    printf("#            tau   Rd  zd zd45 sigU0 sigW0  RsigU  RsigW    rho0"
           "        n0     n0WD    Sigma0  \n");
    printf("#            Gyr   pc  pc   pc  km/s  km/s     pc     pc  Msun/pc^3"
           "  */pc^3   */pc^3  Msun/pc^2\n");
    }
    for (int i = 0; i < 8; i++) {
        int rd    = (i == 0) ? context.density.Rd.data()[0] : (i < 7) ? context.density.Rd.data()[1] : context.density.Rd.data()[2];
        int zdtmp = (context.density.hDISK == 0) ? context.density.zd.data()[i] : context.density.zd45.data()[i];
        double hsigU = (i < 7) ? context.kinematics.hsigUt : context.kinematics.hsigUT;
        double hsigW = (i < 7) ? context.kinematics.hsigWt : context.kinematics.hsigWT;
        double sigW0 = (i < 7) ? context.kinematics.sigW10d * pow((context.kinematics.medtauds.data()[i]+0.01)/10.01, context.kinematics.betaW) : context.kinematics.sigW0td;
        double sigU0 = (i < 7) ? context.kinematics.sigU10d * pow((context.kinematics.medtauds.data()[i]+0.01)/10.01, context.kinematics.betaU) : context.kinematics.sigU0td;
        double SigmaD = 2 * context.density.zd.data()[i] * context.density.rho0d.data()[i];
        if (emit_cli_output)
        printf("#   Disk%d: %5.2f %4d %3.0f  %3d %5.2f %5.2f %6.0f %6.0f"
               "   %.2e %.2e %.2e   %.2e\n",
               i+1, context.kinematics.medtauds.data()[i], rd, context.density.zd.data()[i], zdtmp, sigU0, sigW0, hsigU, hsigW,
               context.density.rho0d.data()[i], context.density.n0d.data()[i], context.density.n0d.data()[i]-context.density.n0MSd.data()[i], SigmaD);
    }
    if (emit_cli_output) {
    printf("#------------ Disk BH scale height --------------\n");
    printf("#            BHhd = %d\n", BHhd);
    if (BHhd == 1)
        printf("#   (Rhd, beta_v) = (%6.1f pc, %6.4f) \n", RhdBH0, betaBH);
    printf("#        UseSigBH = %d\n", UseSigBH);
    printf("#         MXDkick = %d\n", MXDkick);
    printf("#   (vkickNS, vkickBH) = (%6.1f , %6.1f) km/s\n", vkickNS, vkickBH);
    }

    // ------- NSD option parsing -------
    double MND = 0;
    if (fabs(lSIMU) < 5 && fabs(bSIMU) < 2) context.density.ND = 3;
    context.density.ND = getOptiond(argc, argv, "NSD", 1, context.density.ND);
    if (typed_config && typed_config->model.nsd.enabled >= 0) {
        context.density.ND = typed_config->model.nsd.enabled;
    }
    if (context.density.ND > 0) context.density.ND = 3;
    if (context.density.ND == 1) { MND = 2.0e+09; context.density.x0ND = 250; context.density.y0ND = 125; context.density.z0ND = 50; }
    if (context.density.ND == 2) { MND = 7.0e+08; context.density.x0ND =  74; context.density.y0ND =  74; context.density.z0ND = 26; }
    context.density.x0ND = getOptiond(argc, argv, "x0ND", 1, context.density.x0ND);
    context.density.y0ND = getOptiond(argc, argv, "y0ND", 1, context.density.y0ND);
    context.density.z0ND = getOptiond(argc, argv, "z0ND", 1, context.density.z0ND);
    MND  = getOptiond(argc, argv, "MND",  1,  MND);
    if (typed_config) {
        const auto &nsd = typed_config->model.nsd;
        context.density.x0ND = nsd.x0;
        context.density.y0ND = nsd.y0;
        context.density.z0ND = nsd.z0;
        if (nsd.mass > 0.0) MND = nsd.mass;
    }
    if (context.density.ND)
        context.density.rho0ND = (context.density.ND == 3) ? 1 : 0.25*MND/PI/context.density.x0ND/context.density.y0ND/context.density.z0ND;

    // ------- Bulge + disk + NSD normalization -------
    Initializer().finalize_density_normalization(context);
    double massentire = crude_integrate(context, 6000, 3000, 3000, 30) * context.density.rho0b;

    prepared.nsd_moments.initialize_if_enabled(context);
    prepared.nsd_moments_active = true;

    if (emit_cli_output) {
    printf("#--- Bulge: (alpha_bar, Mbar)= (%.1f deg, %.2e Msun) ---\n",
           context.density.thetaD, massentire);
    printf("#   (M_MS, M_REM)ave= (%.6f %.6f) Msun/*, fM_REM= %.4f\n",
           1/context.density.m2nb_MS, 1/context.density.m2nb_WD, 1-context.density.fb_MS);
    printf("#   rho%d: rho0b= %5.2f, (x0,y0,z0,Rc)= (%4.0f,%4.0f,%3.0f,%4.0f),"
           " (C1,C2,C3)= (%.1f,%.1f,%.1f)\n",
           context.density.model, context.density.rho0b, context.density.x0_1, context.density.y0_1, context.density.z0_1, context.density.Rc, context.density.C1, context.density.C2, context.density.C3);
    if (context.density.addX >= 5)
        printf("#     X%d: rho0X= %5.2f, (x0,y0,z0,Rc)= (%4.0f,%4.0f,%3.0f,%4.0f)\n",
               context.density.addX, context.density.rho0b*context.density.fX, context.density.x0_X, context.density.y0_X, context.density.z0_X, context.density.Rc_X);
    printf("#   (Omega_p, vx_str)= (%.1f km/s/kpc, %3.0f[1-e^{-(|yb|/%4.0f)^2}] km/s)\n",
           context.kinematics.Omega_p, context.kinematics.vx_str, context.kinematics.y0_str);
    printf("#   ND= %d  SH= %d\n", context.density.ND, context.density.SH);
    }

    // ------- PA calculation -------
    double PA, cosPA, sinPA;
    void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
    calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);

    // ------- Sampling options -------
    auto &sampling_options = context.sampling;
    Initializer().read_sampling_options(context, argc, argv, cosPA, sinPA);
    if (typed_config) {
        sampling_options.n_simu = typed_config->n_simu;
        sampling_options.n_like_min = typed_config->sampling.n_like_min;
        sampling_options.v_earth_l = typed_config->sampling.v_earth_l;
        sampling_options.v_earth_b = typed_config->sampling.v_earth_b;
        sampling_options.gamma_ds = typed_config->sampling.gamma_ds;
        sampling_options.weight_lens_distance = typed_config->sampling.weight_lens_distance;
        sampling_options.weight_lens_mass = typed_config->sampling.weight_lens_mass;
        sampling_options.no_gamma_importance_sampling = typed_config->sampling.no_gamma_importance_sampling;
        sampling_options.small_gamma = typed_config->sampling.small_gamma;
        sampling_options.verbosity = typed_config->sampling.verbosity;
        sampling_options.uniform_likelihood = typed_config->sampling.uniform_likelihood;
        sampling_options.binary = typed_config->sampling.binary;
        sampling_options.remnant = typed_config->sampling.remnant;
        sampling_options.only_white_dwarf = typed_config->sampling.only_white_dwarf;
        sampling_options.calc_prior_piE = typed_config->sampling.calc_prior_piE;
        sampling_options.calc_prior_thetaE = typed_config->sampling.calc_prior_thetaE;
    }

    // ------- Observation config -------
    prepared.observation = ObservationConfig::from_cli(
        argc, argv, sampling_options.uniform_likelihood);
    if (typed_config) {
        prepared.observation = typed_config->observation;
        if (typed_config->observed_tE != 54.5 ||
            typed_config->observed_tE_error != 99999999999.0) {
            prepared.observation.tE_obs = typed_config->observed_tE;
            prepared.observation.tE_err = typed_config->observed_tE_error;
        }
        const int no_importance_sampling = (int)getOptiond(argc, argv, "NOIS", 1, 0);
        const double default_importance_width =
            (sampling_options.uniform_likelihood == 1) ? 1.02 : 4.0;
        const double tE_importance_width =
            getOptiond(argc, argv, "fIStE", 1, default_importance_width);
        if (no_importance_sampling == 0 &&
            prepared.observation.tE_obs - prepared.observation.tE_err > 0 &&
            prepared.observation.tE_err > 0 &&
            prepared.observation.tE_max - prepared.observation.tE_min == 0) {
            prepared.observation.tE_min = prepared.observation.tE_obs - tE_importance_width * prepared.observation.tE_err;
            prepared.observation.tE_max = prepared.observation.tE_obs + tE_importance_width * prepared.observation.tE_err;
            if (prepared.observation.tE_min <= 0 || prepared.observation.tE_det == 1) prepared.observation.tE_min = 1e-10;
            if (prepared.observation.tE_det == 2) prepared.observation.tE_max = 1e+6;
        }
    }

    // ------- Extinction -------
    if (DMrc == 0)
        DMrc = 14.3955 - 0.0239*lSIMU + 0.0122*fabs(bSIMU) + 0.128;
    if (sampling_options.n_simu == 0) exit(1);

    context.density.lDs    = (double*)malloc(sizeof(double) * 1);
    context.density.bDs    = (double*)malloc(sizeof(double) * 1);
    prepared.density_sightline_allocated = true;
    context.density.lDs[0] = lSIMU;
    context.density.bDs[0] = bSIMU;

    double hdust  = getOptiond(argc, argv, "hdust", 1, 164.0);
    if (typed_config) {
        hdust = typed_config->source.dust_scale_height_pc;
    }
    double cosb   = cos(bSIMU/180.0*PI), sinb = sin(bSIMU/180.0*PI);
    double cosl   = cos(lSIMU/180.0*PI), sinl = sin(lSIMU/180.0*PI);
    double hscale = hdust / (fabs(sinb) + 0.0001);
    double Dmean  = (DMrc > 0) ? pow(10, 0.2*DMrc) * 10 : -9.99;
    gmodel::ExponentialDustExtinction extinction(hscale, Dmean, AIrc, AKrc, EVIrc);
    double AI0  = extinction.ai0();
    double AK0  = extinction.ak0();
    double EVI0 = extinction.evi0();

    if (emit_cli_output) {
    printf("#---------- Input parameters ----------\n");
    printf("#    CenSgrA= %d\n", context.spatial.center_on_sgr_a);
    printf("#    UNIFORM= %d\n", sampling_options.uniform_likelihood);
    printf("#    REMNANT= %d     BINARY= %d\n",
           sampling_options.remnant, sampling_options.binary);
    printf("#   (Nsimu, NlikeMIN)= (%ld, %ld)\n",
           sampling_options.n_simu, sampling_options.n_like_min);
    printf("#          (l, b, PA)= (%6.2f, %6.2f, %5.2f) deg.\n", lSIMU, bSIMU, PA);
    if (prepared.observation.tE_err > 0)
        printf("#     tE = %.3f +- %.3f det= %d\n", prepared.observation.tE_obs, prepared.observation.tE_err, prepared.observation.tE_det);
    if (prepared.observation.thetaE_err > 0)
        printf("# thetaE = %.3f +- %.3f det= %d\n",
               prepared.observation.thetaE_obs, prepared.observation.thetaE_err, prepared.observation.thetaE_det);
    if (prepared.observation.piE_err > 0)
        printf("#    piE = %.3f +- %.3f det= %d\n", prepared.observation.piE_obs, prepared.observation.piE_err, prepared.observation.piE_det);
    if (AI0 > 0)
        printf("#  Consider %.2f < Is < %.2f, (hdust,Dmean,AIrc,AI0)="
               "(%.0f,%.0f,%.2f,%.2f)\n", Isst, Isen, hdust, Dmean, AIrc, AI0);
    if (EVI0 > 0)
        printf("#  Consider %.2f < VIs < %.2f, (hdust,Dmean,EVIrc,EVI0)="
               "(%.0f,%.0f,%.2f,%.2f)\n", VIsst, VIsen, hdust, Dmean, EVIrc, EVI0);
    if (prepared.observation.tE_max - prepared.observation.tE_min > 0)
        printf("#     tErange     : %.4f - %.4f\n", prepared.observation.tE_min, prepared.observation.tE_max);
    if (prepared.observation.thetaE_max - prepared.observation.thetaE_min > 0)
        printf("#     thetaErange : %.4f - %.4f\n", prepared.observation.thetaE_min, prepared.observation.thetaE_max);
    if (prepared.observation.piE_max - prepared.observation.piE_min > 0)
        printf("#     piErange    : %.4f - %.4f\n", prepared.observation.piE_min, prepared.observation.piE_max);
    }

    // ------- wtM_L adjustment -------
    if (sampling_options.weight_lens_mass != 0) {
        alpha0_B += sampling_options.weight_lens_mass;
        alpha1_B += sampling_options.weight_lens_mass;
        alpha2_B += sampling_options.weight_lens_mass;
        alpha3_B += sampling_options.weight_lens_mass;
        alpha4_B += sampling_options.weight_lens_mass;
        store_IMF_nBs(context, 0, logMass_B, PlogM_B, PlogM_cum_norm_B, imptiles_B,
                      M0_B, M1_B, M2_B, M3_B, Ml, Mu,
                      alpha1_B, alpha2_B, alpha3_B, alpha4_B, alpha0_B);
    }

    // ------- Density grid -------
    int Dmax = getOptiond(argc, argv, "Dmax", 1, 16000);
    if (typed_config) {
        Dmax = typed_config->runtime.max_distance_pc;
    }
    int npri = (context.density.ND > 0) ? 40 : 10;
    if (printBHfac) npri = 100000;
    npri = getOptioni(argc, argv, "npri", 1, npri);
    int printrhoS = getOptiond(argc, argv, "printrhoS", 1, 0);
    if (!emit_cli_output) {
        npri = 0;
        printBHfac = 0;
        printrhoS = 0;
    }

    LineOfSightDensityGridConfig grid_cfg;
    grid_cfg.Dmax      = Dmax;
    grid_cfg.l         = lSIMU;
    grid_cfg.b         = bSIMU;
    grid_cfg.BHhd      = BHhd;
    grid_cfg.BHhb      = BHhb;
    grid_cfg.fixRhdBH  = fixRhdBH;
    grid_cfg.vkickBH   = vkickBH;
    grid_cfg.RhdBH0    = RhdBH0;
    grid_cfg.betaBH    = betaBH;
    grid_cfg.UseSigBH  = UseSigBH;
    grid_cfg.AI0       = AI0;
    grid_cfg.EVI0      = EVI0;
    grid_cfg.gammaDs   = sampling_options.gamma_ds;
    grid_cfg.Isst      = Isst;
    grid_cfg.Isen      = Isen;
    grid_cfg.VIsst     = VIsst;
    grid_cfg.VIsen     = VIsen;
    grid_cfg.wtD_L     = sampling_options.weight_lens_distance;
    grid_cfg.check_D   = (CheckD == 1);
    grid_cfg.npri      = npri;
    grid_cfg.printBHfac = (printBHfac != 0);
    grid_cfg.printrhoS = (printrhoS != 0);
    grid_cfg.extinction = &extinction;

    if (emit_cli_output)
        printf("#----- Mass density along (l,b)= (%.3f,%.3f) --------\n",
               lSIMU, bSIMU);
    prepared.grid.build(context, grid_cfg);

    // ------- Optional CALCTAU -------
    int    CALCTAU = getOptioni(argc, argv, "CALCTAU", 1, 0);
    if (typed_config) {
        CALCTAU = typed_config->runtime.calculate_optical_depth;
    }
    double tauall  = 0, Nsall = 0;
    if (CALCTAU && Isen - Isst > 0 && AIrc > 0) {
        void calc_opticaldepth(RunContext &ctx, double *tauall, double *Nsall, int idata,
                               int Dsmax21, double AI0, double hscale,
                               double Isst, double Isen);
        calc_opticaldepth(context, &tauall, &Nsall, 0, Dmax, AI0, hscale, Isst, Isen);
    }

    // ------- CheckD mode -------
    if (CheckD == 1) {
        double getcumu2xist(int n, double *x, double *F, double *f,
                            double Freq, int ist, int inv);
        for (int i = 0; i < 500000; i++) {
            double ran = context.runtime.rng->uniform(), cumu = 0;
            int j_L, j_S;
            for (j_L = 0; j_L < context.density.ncomp; j_L++) {
                cumu += prepared.grid.cumu_rho_L_data(j_L)[prepared.grid.nbin()] /
                        prepared.grid.cumu_rho_all_L_data()[prepared.grid.nbin()];
                if (ran < cumu) break;
            }
            cumu = 0;
            for (j_S = 0; j_S < context.density.ncomp; j_S++) {
                cumu += prepared.grid.cumu_rho_S_data(j_S)[prepared.grid.nbin()] /
                        prepared.grid.cumu_rho_all_S_data()[prepared.grid.nbin()];
                if (ran < cumu) break;
            }
            double d_L = prepared.grid.sample_lens_distance(j_L, 0, prepared.grid.nbin(), context.runtime.rng->uniform());
            double d_S = prepared.grid.sample_source_distance(j_S, context.runtime.rng->uniform());
            void get_vxyz_ran(RunContext &ctx, double *vxyz, int i, double tau, double D, double lD, double bD);
            double vxyz_L[3] = {};
            double tau_l = (j_L == 9) ? context.stellar.mageND + context.stellar.sageND*context.runtime.rng->gaussian()
                         : (j_L == 8) ? context.stellar.mageB + context.stellar.sageB*context.runtime.rng->gaussian()
                         : (j_L == 10) ? 14.0 : context.kinematics.medtauds.data()[j_L];
            get_vxyz_ran(context, vxyz_L, j_L, tau_l, d_L, context.density.lDs[0], context.density.bDs[0]);
            double v_l = sqrt(vxyz_L[0]*vxyz_L[0] + vxyz_L[1]*vxyz_L[1]
                            + vxyz_L[2]*vxyz_L[2]);
            if (emit_cli_output)
                printf("%d %6.0f %d %6.0f %.6f %7.2f %7.2f %7.2f %7.2f\n",
                   j_S, d_S, j_L, d_L, sqrt(d_L/8000.0),
                   vxyz_L[0], vxyz_L[1], vxyz_L[2], v_l);
        }
        exit(1);
    }

    // Release luminosity tables (no longer needed after grid build)
    prepared.population.release_luminosity_functions(context);
    prepared.luminosity_functions_active = false;

    // ------- MassFunction wrapper -------
    prepared.mass_function.init_from_population(prepared.population, context);

    // ------- Monte Carlo simulation -------
    if (emit_cli_output)
        printf("#----- Output of Monte Carlo simulation w/ VERBOSITY= %d and seed= %ld ----\n",
               sampling_options.verbosity, context.seed);

    prepared.event_config.cosPA   = cosPA;    prepared.event_config.sinPA = sinPA;
    prepared.event_config.cosb    = cosb;     prepared.event_config.sinb  = sinb;
    prepared.event_config.cosl    = cosl;     prepared.event_config.sinl  = sinl;
    prepared.event_config.l       = lSIMU;   prepared.event_config.b     = bSIMU;
    prepared.event_config.extinction = &extinction;
    prepared.event_config.AI0     = AI0;     prepared.event_config.AK0  = AK0;
    prepared.event_config.BHhd    = BHhd;    prepared.event_config.BHhb = BHhb;
    prepared.event_config.vkickBH = vkickBH; prepared.event_config.vkickNS = vkickNS;
    prepared.event_config.MXDkick = MXDkick; prepared.event_config.betaBH  = betaBH;
    prepared.event_config.nallS   = prepared.grid.total_source_count();
    prepared.event_config.Nsall   = Nsall;
    prepared.event_config.tauall  = tauall;
    prepared.event_config.Dmean   = Dmean;
    if (typed_config && typed_config->forward_source.enabled != 0) {
        prepared.forward_source_generator =
            std::make_unique<gmodel::ForwardSourceGenerator>(make_forward_source_generator(*typed_config));
        prepared.forward_source_rng =
            std::make_unique<RandomEngine>(typed_config->seed + 0x9e3779b9UL);
        prepared.event_config.forward_source_generator = prepared.forward_source_generator.get();
        prepared.event_config.forward_source_rng = prepared.forward_source_rng.get();
        prepared.event_config.source_min_initial_mass_msun =
            typed_config->forward_source.min_initial_mass_msun;
        prepared.event_config.source_max_initial_mass_msun =
            typed_config->forward_source.max_initial_mass_msun;
        prepared.event_config.match_source_selection =
            typed_config->forward_source.match_source_selection != 0;
        prepared.event_config.max_source_selection_attempts =
            typed_config->forward_source.max_selection_attempts;
        prepared.event_config.source_i_min = Isst;
        prepared.event_config.source_i_max = Isen;
        prepared.event_config.source_vi_min = VIsst;
        prepared.event_config.source_vi_max = VIsen;
        prepared.event_config.source_selection_bands =
            typed_config->forward_source.selection_bands;
        prepared.event_config.source_selection_min_magnitudes =
            typed_config->forward_source.selection_min_magnitudes;
        prepared.event_config.source_selection_max_magnitudes =
            typed_config->forward_source.selection_max_magnitudes;
        prepared.event_config.source_selection_apparent_magnitudes =
            typed_config->forward_source.selection_apparent_magnitudes;
    }

    const int code = run_prepared_events(context, prepared, std::move(custom_likelihood),
                                         std::move(event_sink), emit_cli_output);

    // ------- Cleanup -------
    prepared.cleanup(context);
    return code;
}

int Sampler::run_cli(RunContext &context, int argc, char **argv)
{
    return run_sampler_impl(context, nullptr, argc, argv, {}, {}, true);
}

SimulationResult Sampler::simulate(RunContext &context, int argc, char **argv,
                                   LikelihoodFunction likelihood)
{
    SimulationResult result;
    const int code = run_sampler_impl(
        context, nullptr, argc, argv, std::move(likelihood),
        [&result](const Event &event) {
            result.events.push_back(event);
        },
        false);
    if (code != 0) {
        throw std::runtime_error("genulens scientific backend returned non-zero status");
    }
    return result;
}

SimulationResult Sampler::simulate(RunContext &context, const GenulensConfig &config,
                                   int argc, char **argv,
                                   LikelihoodFunction likelihood)
{
    SimulationResult result;
    result.verbosity = config.sampling.verbosity;
    result.include_source_properties = config.forward_source.enabled != 0;
    result.source_property_bands = forward_source_bands_for_config(config);
    const int code = run_sampler_impl(
        context, &config, argc, argv, std::move(likelihood),
        [&result](const Event &event) {
            result.events.push_back(event);
        },
        false);
    if (code != 0) {
        throw std::runtime_error("genulens scientific backend returned non-zero status");
    }
    return result;
}

} // namespace genulens
