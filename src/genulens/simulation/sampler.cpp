/* Generate microlensing events following the Galactic model developed by Koshimoto, Baba & Bennett (2021).
 * N. Koshimoto wrote the original .c version and C. Ranc converted it into .cpp.
 * This is version 1.2.1+ of genulens (refactored). */
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "genulens/cli/option.h"
#include "genulens/io/input_data.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/extinction.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/rng.hpp"
#include "genulens/simulation/event_sampler.hpp"
#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/los_density_grid.hpp"
#include "genulens/simulation/mass_function.hpp"
#include "genulens/simulation/observation_config.hpp"
#include "genulens/simulation/sampler.hpp"
#include "genulens/simulation/velocity_distribution.hpp"

#define fopen(path, mode) genulens::open_input_file((path), (mode))

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"

#define PI 3.1415926535897932385

namespace genulens {




int Sampler::run_cli(RunContext &context, int argc, char **argv)
{
    char curdir[200];
    getcwd(curdir, 200);

    // ------- Model options -------
    int CheckD = getOptiond(argc, argv, "CheckD", 1, 0);
    Initializer().read_model_options(context, argc, argv);
    double M0_B = context.imf_options.m0, M1_B = context.imf_options.m1;
    double M2_B = context.imf_options.m2, M3_B = context.imf_options.m3;
    double Ml   = context.imf_options.ml, Mu   = context.imf_options.mu;
    double alpha0_B = context.imf_options.alpha0, alpha1_B = context.imf_options.alpha1;
    double alpha2_B = context.imf_options.alpha2, alpha3_B = context.imf_options.alpha3;
    double alpha4_B = context.imf_options.alpha4;
    Initializer().finalize_spatial_model(context, argc, argv);

    // ------- Sightline & source-selection -------
    double lSIMU = getOptiond(argc, argv, "l",        1,  1.0);
    double bSIMU = getOptiond(argc, argv, "b",        1, -3.9);
    double Isst  = getOptiond(argc, argv, "Isrange",  1, 14.0);
    double Isen  = getOptiond(argc, argv, "Isrange",  2, 21.0);
    double VIsst = getOptiond(argc, argv, "VIsrange", 1,  0.0);
    double VIsen = getOptiond(argc, argv, "VIsrange", 2,  0.0);
    double AIrc  = getOptiond(argc, argv, "AIrc",     1,  0);
    double EVIrc = getOptiond(argc, argv, "EVIrc",    1,  0);
    double DMrc  = getOptiond(argc, argv, "DMrc",     1,  0);
    double AKrc  = getOptiond(argc, argv, "AKrc",     1,  0);
    if (fabs(lSIMU) > 10 || fabs(bSIMU) > 7 || fabs(bSIMU) < 1.5)
        printf("# WARNING: genulens is designed for |l| < ~10 and ~1.5 < |b| < ~7"
               " and (l,b)= (%.3f, %.3f) is outside the range.\n", lSIMU, bSIMU);

    // ------- Population runtime -------
    PopulationRuntime population;
    population.initialize_mass_function(context, context.imf_options);
    population.initialize_luminosity_functions(context, Isst, Isen, VIsst, VIsen, AIrc, EVIrc);
    population.read_empirical_mass_luminosity();
    double *logMass_B        = population.log_mass;
    double *PlogM_B          = population.mass_probability;
    double *PlogM_cum_norm_B = population.mass_cumulative;
    int    *imptiles_B       = population.mass_percentiles;

    // ------- Kinematic tables -------
    KinematicRuntimeTables kinematic_tables;
    kinematic_tables.initialize_shu_distribution(context);

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

    // ------- Header printing -------
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
    for (int i = 0; i < 8; i++) {
        int rd    = (i == 0) ? context.density.Rd.data()[0] : (i < 7) ? context.density.Rd.data()[1] : context.density.Rd.data()[2];
        int zdtmp = (context.density.hDISK == 0) ? context.density.zd.data()[i] : context.density.zd45.data()[i];
        double hsigU = (i < 7) ? context.kinematics.hsigUt : context.kinematics.hsigUT;
        double hsigW = (i < 7) ? context.kinematics.hsigWt : context.kinematics.hsigWT;
        double sigW0 = (i < 7) ? context.kinematics.sigW10d * pow((context.kinematics.medtauds.data()[i]+0.01)/10.01, context.kinematics.betaW) : context.kinematics.sigW0td;
        double sigU0 = (i < 7) ? context.kinematics.sigU10d * pow((context.kinematics.medtauds.data()[i]+0.01)/10.01, context.kinematics.betaU) : context.kinematics.sigU0td;
        double SigmaD = 2 * context.density.zd.data()[i] * context.density.rho0d.data()[i];
        printf("#   Disk%d: %5.2f %4d %3.0f  %3d %5.2f %5.2f %6.0f %6.0f"
               "   %.2e %.2e %.2e   %.2e\n",
               i+1, context.kinematics.medtauds.data()[i], rd, context.density.zd.data()[i], zdtmp, sigU0, sigW0, hsigU, hsigW,
               context.density.rho0d.data()[i], context.density.n0d.data()[i], context.density.n0d.data()[i]-context.density.n0MSd.data()[i], SigmaD);
    }
    printf("#------------ Disk BH scale height --------------\n");
    printf("#            BHhd = %d\n", BHhd);
    if (BHhd == 1)
        printf("#   (Rhd, beta_v) = (%6.1f pc, %6.4f) \n", RhdBH0, betaBH);
    printf("#        UseSigBH = %d\n", UseSigBH);
    printf("#         MXDkick = %d\n", MXDkick);
    printf("#   (vkickNS, vkickBH) = (%6.1f , %6.1f) km/s\n", vkickNS, vkickBH);

    // ------- NSD option parsing -------
    double MND = 0;
    if (fabs(lSIMU) < 5 && fabs(bSIMU) < 2) context.density.ND = 3;
    context.density.ND = getOptiond(argc, argv, "NSD", 1, context.density.ND);
    if (context.density.ND > 0) context.density.ND = 3;
    if (context.density.ND == 1) { MND = 2.0e+09; context.density.x0ND = 250; context.density.y0ND = 125; context.density.z0ND = 50; }
    if (context.density.ND == 2) { MND = 7.0e+08; context.density.x0ND =  74; context.density.y0ND =  74; context.density.z0ND = 26; }
    context.density.x0ND = getOptiond(argc, argv, "x0ND", 1, context.density.x0ND);
    context.density.y0ND = getOptiond(argc, argv, "y0ND", 1, context.density.y0ND);
    context.density.z0ND = getOptiond(argc, argv, "z0ND", 1, context.density.z0ND);
    MND  = getOptiond(argc, argv, "MND",  1,  MND);
    if (context.density.ND)
        context.density.rho0ND = (context.density.ND == 3) ? 1 : 0.25*MND/PI/context.density.x0ND/context.density.y0ND/context.density.z0ND;

    // ------- Bulge + disk + NSD normalization -------
    Initializer().finalize_density_normalization(context);
    double massentire = crude_integrate(context, 6000, 3000, 3000, 30) * context.density.rho0b;

    NsdMomentRuntime nsd_moments;
    nsd_moments.initialize_if_enabled(context);

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

    // ------- PA calculation -------
    double PA, cosPA, sinPA;
    void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
    calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);

    // ------- Sampling options -------
    auto &sampling_options = context.sampling;
    Initializer().read_sampling_options(context, argc, argv, cosPA, sinPA);

    // ------- Observation config -------
    ObservationConfig obs = ObservationConfig::from_cli(
        argc, argv, sampling_options.uniform_likelihood);

    // ------- Extinction -------
    if (DMrc == 0)
        DMrc = 14.3955 - 0.0239*lSIMU + 0.0122*fabs(bSIMU) + 0.128;
    if (sampling_options.n_simu == 0) exit(1);

    context.density.lDs    = (double*)malloc(sizeof(double) * 1);
    context.density.bDs    = (double*)malloc(sizeof(double) * 1);
    context.density.lDs[0] = lSIMU;
    context.density.bDs[0] = bSIMU;

    double hdust  = getOptiond(argc, argv, "hdust", 1, 164.0);
    double cosb   = cos(bSIMU/180.0*PI), sinb = sin(bSIMU/180.0*PI);
    double cosl   = cos(lSIMU/180.0*PI), sinl = sin(lSIMU/180.0*PI);
    double hscale = hdust / (fabs(sinb) + 0.0001);
    double Dmean  = (DMrc > 0) ? pow(10, 0.2*DMrc) * 10 : -9.99;
    gmodel::ExponentialDustExtinction extinction(hscale, Dmean, AIrc, AKrc, EVIrc);
    double AI0  = extinction.ai0();
    double AK0  = extinction.ak0();
    double EVI0 = extinction.evi0();

    printf("#---------- Input parameters ----------\n");
    printf("#    CenSgrA= %d\n", context.spatial.center_on_sgr_a);
    printf("#    UNIFORM= %d\n", sampling_options.uniform_likelihood);
    printf("#    REMNANT= %d     BINARY= %d\n",
           sampling_options.remnant, sampling_options.binary);
    printf("#   (Nsimu, NlikeMIN)= (%ld, %ld)\n",
           sampling_options.n_simu, sampling_options.n_like_min);
    printf("#          (l, b, PA)= (%6.2f, %6.2f, %5.2f) deg.\n", lSIMU, bSIMU, PA);
    if (obs.tE_err > 0)
        printf("#     tE = %.3f +- %.3f det= %d\n", obs.tE_obs, obs.tE_err, obs.tE_det);
    if (obs.thetaE_err > 0)
        printf("# thetaE = %.3f +- %.3f det= %d\n",
               obs.thetaE_obs, obs.thetaE_err, obs.thetaE_det);
    if (obs.piE_err > 0)
        printf("#    piE = %.3f +- %.3f det= %d\n", obs.piE_obs, obs.piE_err, obs.piE_det);
    if (AI0 > 0)
        printf("#  Consider %.2f < Is < %.2f, (hdust,Dmean,AIrc,AI0)="
               "(%.0f,%.0f,%.2f,%.2f)\n", Isst, Isen, hdust, Dmean, AIrc, AI0);
    if (EVI0 > 0)
        printf("#  Consider %.2f < VIs < %.2f, (hdust,Dmean,EVIrc,EVI0)="
               "(%.0f,%.0f,%.2f,%.2f)\n", VIsst, VIsen, hdust, Dmean, EVIrc, EVI0);
    if (obs.tE_max - obs.tE_min > 0)
        printf("#     tErange     : %.4f - %.4f\n", obs.tE_min, obs.tE_max);
    if (obs.thetaE_max - obs.thetaE_min > 0)
        printf("#     thetaErange : %.4f - %.4f\n", obs.thetaE_min, obs.thetaE_max);
    if (obs.piE_max - obs.piE_min > 0)
        printf("#     piErange    : %.4f - %.4f\n", obs.piE_min, obs.piE_max);

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
    int npri = (context.density.ND > 0) ? 40 : 10;
    if (printBHfac) npri = 100000;
    npri = getOptioni(argc, argv, "npri", 1, npri);
    int printrhoS = getOptiond(argc, argv, "printrhoS", 1, 0);

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

    printf("#----- Mass density along (l,b)= (%.3f,%.3f) --------\n",
           lSIMU, bSIMU);
    LineOfSightDensityGrid grid;
    grid.build(context, grid_cfg);

    // ------- Optional CALCTAU -------
    int    CALCTAU = getOptioni(argc, argv, "CALCTAU", 1, 0);
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
                cumu += grid.cumu_rho_L_data(j_L)[grid.nbin()] /
                        grid.cumu_rho_all_L_data()[grid.nbin()];
                if (ran < cumu) break;
            }
            cumu = 0;
            for (j_S = 0; j_S < context.density.ncomp; j_S++) {
                cumu += grid.cumu_rho_S_data(j_S)[grid.nbin()] /
                        grid.cumu_rho_all_S_data()[grid.nbin()];
                if (ran < cumu) break;
            }
            double d_L = grid.sample_lens_distance(j_L, 0, grid.nbin(), context.runtime.rng->uniform());
            double d_S = grid.sample_source_distance(j_S, context.runtime.rng->uniform());
            void get_vxyz_ran(RunContext &ctx, double *vxyz, int i, double tau, double D, double lD, double bD);
            double vxyz_L[3] = {};
            double tau_l = (j_L == 9) ? context.stellar.mageND + context.stellar.sageND*context.runtime.rng->gaussian()
                         : (j_L == 8) ? context.stellar.mageB + context.stellar.sageB*context.runtime.rng->gaussian()
                         : (j_L == 10) ? 14.0 : context.kinematics.medtauds.data()[j_L];
            get_vxyz_ran(context, vxyz_L, j_L, tau_l, d_L, context.density.lDs[0], context.density.bDs[0]);
            double v_l = sqrt(vxyz_L[0]*vxyz_L[0] + vxyz_L[1]*vxyz_L[1]
                            + vxyz_L[2]*vxyz_L[2]);
            printf("%d %6.0f %d %6.0f %.6f %7.2f %7.2f %7.2f %7.2f\n",
                   j_S, d_S, j_L, d_L, sqrt(d_L/8000.0),
                   vxyz_L[0], vxyz_L[1], vxyz_L[2], v_l);
        }
        exit(1);
    }

    // Release luminosity tables (no longer needed after grid build)
    population.release_luminosity_functions(context);

    // ------- MassFunction wrapper -------
    MassFunction mf;
    mf.init_from_population(population, context);

    // ------- Monte Carlo simulation -------
    printf("#----- Output of Monte Carlo simulation w/ VERBOSITY= %d and seed= %ld ----\n",
           sampling_options.verbosity, context.seed);

    EventSampler::Config es_cfg;
    es_cfg.cosPA   = cosPA;    es_cfg.sinPA = sinPA;
    es_cfg.cosb    = cosb;     es_cfg.sinb  = sinb;
    es_cfg.cosl    = cosl;     es_cfg.sinl  = sinl;
    es_cfg.l       = lSIMU;   es_cfg.b     = bSIMU;
    es_cfg.extinction = &extinction;
    es_cfg.AI0     = AI0;     es_cfg.AK0  = AK0;
    es_cfg.BHhd    = BHhd;    es_cfg.BHhb = BHhb;
    es_cfg.vkickBH = vkickBH; es_cfg.vkickNS = vkickNS;
    es_cfg.MXDkick = MXDkick; es_cfg.betaBH  = betaBH;
    es_cfg.nallS   = grid.total_source_count();
    es_cfg.Nsall   = Nsall;
    es_cfg.tauall  = tauall;
    es_cfg.Dmean   = Dmean;

    EventSampler event_sampler;
    event_sampler.run_cli(context, grid, population, mf, es_cfg, obs);

    // ------- Cleanup -------
    nsd_moments.release_if_enabled(context);
    population.release_all(context);
    free(context.density.lDs);
    free(context.density.bDs);
    kinematic_tables.release_all(context);
    return 0;
}

} // namespace genulens
