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

genulens::RunContext *active_state = nullptr;

double ran1()   { return active_state->runtime.rng->uniform(); }
double gasdev() { return active_state->runtime.rng->gaussian(); }

// All active_state macros (unchanged from original, still needed by runtime helpers)
#define seed active_state->seed
#define B14disk active_state->B14disk
#define B14vbar active_state->B14vbar
#define xyzSgrA active_state->xyzSgrA.data()
#define ncomp active_state->density.ncomp
#define tSFR active_state->stellar.tSFR
#define rhot0 active_state->density.rhot0
#define MiniWDmax active_state->stellar.MiniWDmax
#define agesD active_state->stellar.agesD.data()
#define agesB active_state->stellar.agesB.data()
#define agesND active_state->stellar.agesND.data()
#define MinidieD active_state->stellar.MinidieD.data()
#define MinidieB active_state->stellar.MinidieB.data()
#define MinidieND active_state->stellar.MinidieND.data()
#define nageD active_state->stellar.nageD
#define nageB active_state->stellar.nageB
#define nageND active_state->stellar.nageND
#define mageB active_state->stellar.mageB
#define sageB active_state->stellar.sageB
#define mageND active_state->stellar.mageND
#define sageND active_state->stellar.sageND
#define nm active_state->stellar.nm
#define logMst active_state->stellar.logMst
#define dlogM active_state->stellar.dlogM
#define fb_MS active_state->density.fb_MS
#define m2nb_MS active_state->density.m2nb_MS
#define m2nb_WD active_state->density.m2nb_WD
#define nMS2nRGb active_state->density.nMS2nRGb
#define rho0b active_state->density.rho0b
#define n0MSb active_state->density.n0MSb
#define n0RGb active_state->density.n0RGb
#define n0b active_state->density.n0b
#define ND active_state->density.ND
#define x0ND active_state->density.x0ND
#define y0ND active_state->density.y0ND
#define z0ND active_state->density.z0ND
#define C1ND active_state->density.C1ND
#define rho0ND active_state->density.rho0ND
#define n0MSND active_state->density.n0MSND
#define n0RGND active_state->density.n0RGND
#define n0ND active_state->density.n0ND
#define fND_MS active_state->density.fND_MS
#define m2nND_MS active_state->density.m2nND_MS
#define m2nND_WD active_state->density.m2nND_WD
#define nMS2nRGND active_state->density.nMS2nRGND
#define SH active_state->density.SH
#define rho0SHMS active_state->density.rho0SHMS
#define epsSH active_state->density.epsSH
#define alphaSH active_state->density.alphaSH
#define acSH2 active_state->density.acSH2
#define rho0SH active_state->density.rho0SH
#define n0MSSH active_state->density.n0MSSH
#define n0RGSH active_state->density.n0RGSH
#define n0SH active_state->density.n0SH
#define rho0d active_state->density.rho0d.data()
#define n0d active_state->density.n0d.data()
#define n0MSd active_state->density.n0MSd.data()
#define n0RGd active_state->density.n0RGd.data()
#define y0d active_state->density.y0d.data()
#define Rd active_state->density.Rd.data()
#define Rh active_state->density.Rh
#define Rdbreak active_state->density.Rdbreak
#define nh active_state->density.nh
#define zd active_state->density.zd.data()
#define zd45 active_state->density.zd45.data()
#define DISK active_state->density.DISK
#define hDISK active_state->density.hDISK
#define addX active_state->density.addX
#define model active_state->density.model
#define R0 active_state->density.R0
#define thetaD active_state->density.thetaD
#define x0_1 active_state->density.x0_1
#define y0_1 active_state->density.y0_1
#define z0_1 active_state->density.z0_1
#define C1 active_state->density.C1
#define C2 active_state->density.C2
#define C3 active_state->density.C3
#define Rc active_state->density.Rc
#define frho0b active_state->density.frho0b
#define costheta active_state->density.costheta
#define sintheta active_state->density.sintheta
#define zb_c active_state->density.zb_c
#define x0_X active_state->density.x0_X
#define y0_X active_state->density.y0_X
#define z0_X active_state->density.z0_X
#define C1_X active_state->density.C1_X
#define C2_X active_state->density.C2_X
#define b_zX active_state->density.b_zX
#define fX active_state->density.fX
#define Rsin active_state->density.Rsin
#define b_zY active_state->density.b_zY
#define Rc_X active_state->density.Rc_X
#define lDs active_state->density.lDs
#define bDs active_state->density.bDs
#define nMIs active_state->luminosity.nMIs
#define nVIs active_state->luminosity.nVIs
#define MIs active_state->luminosity.MIs
#define CumuN_MIs active_state->luminosity.CumuN_MIs
#define dILF active_state->luminosity.dILF
#define VIs active_state->luminosity.VIs
#define f_VI_Is active_state->luminosity.f_VI_Is
#define dVILF active_state->luminosity.dVILF
#define nVcs active_state->kinematics.nVcs
#define Rcs active_state->kinematics.Rcs.data()
#define Vcs active_state->kinematics.Vcs.data()
#define vxsun active_state->kinematics.vxsun
#define Vsun active_state->kinematics.Vsun
#define vzsun active_state->kinematics.vzsun
#define vysun active_state->kinematics.vysun
#define fgsShu active_state->kinematics.fgsShu
#define PRRgShus active_state->kinematics.PRRgShus
#define cumu_PRRgs active_state->kinematics.cumu_PRRgs
#define n_fgsShu active_state->kinematics.n_fgsShu
#define kptiles active_state->kinematics.kptiles
#define hsigUt active_state->kinematics.hsigUt
#define hsigWt active_state->kinematics.hsigWt
#define hsigUT active_state->kinematics.hsigUT
#define hsigWT active_state->kinematics.hsigWT
#define betaU active_state->kinematics.betaU
#define betaW active_state->kinematics.betaW
#define sigU10d active_state->kinematics.sigU10d
#define sigW10d active_state->kinematics.sigW10d
#define sigU0td active_state->kinematics.sigU0td
#define sigW0td active_state->kinematics.sigW0td
#define medtauds active_state->kinematics.medtauds.data()
#define zstShu active_state->kinematics.zstShu
#define zenShu active_state->kinematics.zenShu
#define dzShu active_state->kinematics.dzShu
#define RstShu active_state->kinematics.RstShu
#define RenShu active_state->kinematics.RenShu
#define dRShu active_state->kinematics.dRShu
#define model_vb active_state->kinematics.model_vb
#define model_vbz active_state->kinematics.model_vbz
#define Omega_p active_state->kinematics.Omega_p
#define x0_vb active_state->kinematics.x0_vb
#define y0_vb active_state->kinematics.y0_vb
#define z0_vb active_state->kinematics.z0_vb
#define C1_vb active_state->kinematics.C1_vb
#define C2_vb active_state->kinematics.C2_vb
#define C3_vb active_state->kinematics.C3_vb
#define sigx_vb active_state->kinematics.sigx_vb
#define sigy_vb active_state->kinematics.sigy_vb
#define sigz_vb active_state->kinematics.sigz_vb
#define vx_str active_state->kinematics.vx_str
#define y0_str active_state->kinematics.y0_str
#define sigx_vb0 active_state->kinematics.sigx_vb0
#define sigy_vb0 active_state->kinematics.sigy_vb0
#define sigz_vb0 active_state->kinematics.sigz_vb0
#define x0_vbz active_state->kinematics.x0_vbz
#define y0_vbz active_state->kinematics.y0_vbz
#define z0_vbz active_state->kinematics.z0_vbz
#define C1_vbz active_state->kinematics.C1_vbz
#define C2_vbz active_state->kinematics.C2_vbz
#define C3_vbz active_state->kinematics.C3_vbz
#define sigU_SH active_state->kinematics.sigU_SH
#define sigV_SH active_state->kinematics.sigV_SH
#define sigW_SH active_state->kinematics.sigW_SH
#define logrhoNDs active_state->nsd_moments.logrhoNDs
#define vphiNDs active_state->nsd_moments.vphiNDs
#define logsigvNDs active_state->nsd_moments.logsigvNDs
#define corRzNDs active_state->nsd_moments.corRzNDs
#define zstND active_state->nsd_moments.zstND
#define zenND active_state->nsd_moments.zenND
#define dzND active_state->nsd_moments.dzND
#define RstND active_state->nsd_moments.RstND
#define RenND active_state->nsd_moments.RenND
#define dRND active_state->nsd_moments.dRND
#define nzND active_state->nsd_moments.nzND
#define nRND active_state->nsd_moments.nRND

int Sampler::run_cli(RunContext &context, int argc, char **argv)
{
    active_state = &context;
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
    population.initialize_mass_function(context.imf_options);
    population.initialize_luminosity_functions(Isst, Isen, VIsst, VIsen, AIrc, EVIrc);
    population.read_empirical_mass_luminosity();
    double *logMass_B        = population.log_mass;
    double *PlogM_B          = population.mass_probability;
    double *PlogM_cum_norm_B = population.mass_cumulative;
    int    *imptiles_B       = population.mass_percentiles;

    // ------- Kinematic tables -------
    KinematicRuntimeTables kinematic_tables;
    kinematic_tables.initialize_shu_distribution();

    VelocityDistribution vel_dist;
    vel_dist.bind(context);

    // ------- Disk normalization -------
    y0d[0] = (DISK == 1) ? exp(-R0/Rd[0] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[0]);
    y0d[1] = (DISK == 1) ? exp(-R0/Rd[1] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[1]);
    y0d[2] = (DISK == 1) ? exp(-R0/Rd[2] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[2]);

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
    printf("#         (R, z)sun= (%5.0f, %5.0f) pc\n", R0, 25.0);
    printf("#   (vx, vy, vz)sun= (%5.1f, %5.1f, %4.1f) km/s\n", vxsun, vysun, vzsun);
    printf("#------------ Disk model: (DISK, hDISK, tSFR)= ( %d , %d , %.1f Gyr ) -----\n",
           DISK, hDISK, tSFR);
    printf("#            tau   Rd  zd zd45 sigU0 sigW0  RsigU  RsigW    rho0"
           "        n0     n0WD    Sigma0  \n");
    printf("#            Gyr   pc  pc   pc  km/s  km/s     pc     pc  Msun/pc^3"
           "  */pc^3   */pc^3  Msun/pc^2\n");
    double MVVVd = 0, Mind = 0;
    for (int i = 0; i < 8; i++) {
        int rd    = (i == 0) ? Rd[0] : (i < 7) ? Rd[1] : Rd[2];
        int zdtmp = (hDISK == 0) ? zd[i] : zd45[i];
        double MVVVtmp = rho0d[i] * exp((R0 - Rdbreak)/rd) * 2200*2 * 1400*2 * zd[i] / zdtmp;
        double ztmp = 1200.0/zdtmp;
        Mind    += 2 * zdtmp * MVVVtmp / 4400 / 2800 * PI * Rdbreak * Rdbreak;
        MVVVtmp *= (i < 7) ? 2 * zdtmp * (exp(2*ztmp)-1)/(exp(2*ztmp)+1)
                           : 2 * zdtmp * (1 - exp(-ztmp));
        MVVVd += MVVVtmp;
        double hsigU = (i < 7) ? hsigUt : hsigUT;
        double hsigW = (i < 7) ? hsigWt : hsigWT;
        double sigW0 = (i < 7) ? sigW10d * pow((medtauds[i]+0.01)/10.01, betaW) : sigW0td;
        double sigU0 = (i < 7) ? sigU10d * pow((medtauds[i]+0.01)/10.01, betaU) : sigU0td;
        double SigmaD = 2 * zd[i] * rho0d[i];
        printf("#   Disk%d: %5.2f %4d %3.0f  %3d %5.2f %5.2f %6.0f %6.0f"
               "   %.2e %.2e %.2e   %.2e\n",
               i+1, medtauds[i], rd, zd[i], zdtmp, sigU0, sigW0, hsigU, hsigW,
               rho0d[i], n0d[i], n0d[i]-n0MSd[i], SigmaD);
    }
    printf("#------------ Disk BH scale height --------------\n");
    printf("#            BHhd = %d\n", BHhd);
    if (BHhd == 1)
        printf("#   (Rhd, beta_v) = (%6.1f pc, %6.4f) \n", RhdBH0, betaBH);
    printf("#        UseSigBH = %d\n", UseSigBH);
    printf("#         MXDkick = %d\n", MXDkick);
    printf("#   (vkickNS, vkickBH) = (%6.1f , %6.1f) km/s\n", vkickNS, vkickBH);

    // ------- Bulge normalization -------
    double crude_integrate(double xmax, double ymax, double zmax, int nbun);
    double massVVVbox = crude_integrate(2200, 1400, 1200, 15);
    double massentire = crude_integrate(6000, 3000, 3000, 30);
    double fm1 = 1, fmX = 0;
    if (addX >= 5) {
        int addXtmp = addX; addX = 0;
        double mass1all = crude_integrate(6000, 3000, 3000, 30);
        addX = addXtmp;
        fm1 = mass1all / massentire;
        fmX = 1 - fm1;
    }
    double MVVVP17 = 1.32e+10;
    rho0b  = (frho0b * MVVVP17 - MVVVd) / massVVVbox;
    n0MSb  = rho0b * fb_MS * m2nb_MS;
    n0RGb  = n0MSb * nMS2nRGb;
    n0b    = n0MSb + rho0b * (1 - fb_MS) * m2nb_WD;
    massVVVbox *= rho0b;
    massentire *= rho0b;

    // ------- NSD -------
    double MND = 0;
    if (fabs(lSIMU) < 5 && fabs(bSIMU) < 2) ND = 3;
    ND = getOptiond(argc, argv, "NSD", 1, ND);
    if (ND > 0) ND = 3;
    if (ND == 1) { MND = 2.0e+09; x0ND = 250; y0ND = 125; z0ND = 50; }
    if (ND == 2) { MND = 7.0e+08; x0ND =  74; y0ND =  74; z0ND = 26; }
    x0ND = getOptiond(argc, argv, "x0ND", 1, x0ND);
    y0ND = getOptiond(argc, argv, "y0ND", 1, y0ND);
    z0ND = getOptiond(argc, argv, "z0ND", 1, z0ND);
    MND  = getOptiond(argc, argv, "MND",  1,  MND);
    if (ND) {
        rho0ND = (ND == 3) ? 1 : 0.25*MND/PI/x0ND/y0ND/z0ND;
        n0MSND = rho0ND * fND_MS * m2nND_MS;
        n0RGND = n0MSND * nMS2nRGND;
        n0ND   = n0MSND + rho0ND * (1 - fND_MS) * m2nND_WD;
    }
    NsdMomentRuntime nsd_moments;
    nsd_moments.initialize_if_enabled();

    printf("#--- Bulge: (alpha_bar, Mbar)= (%.1f deg, %.2e Msun) ---\n",
           thetaD, massentire);
    printf("#   (M_MS, M_REM)ave= (%.6f %.6f) Msun/*, fM_REM= %.4f\n",
           1/m2nb_MS, 1/m2nb_WD, 1-fb_MS);
    printf("#   rho%d: rho0b= %5.2f, (x0,y0,z0,Rc)= (%4.0f,%4.0f,%3.0f,%4.0f),"
           " (C1,C2,C3)= (%.1f,%.1f,%.1f)\n",
           model, rho0b, x0_1, y0_1, z0_1, Rc, C1, C2, C3);
    if (addX >= 5)
        printf("#     X%d: rho0X= %5.2f, (x0,y0,z0,Rc)= (%4.0f,%4.0f,%3.0f,%4.0f)\n",
               addX, rho0b*fX, x0_X, y0_X, z0_X, Rc_X);
    printf("#   (Omega_p, vx_str)= (%.1f km/s/kpc, %3.0f[1-e^{-(|yb|/%4.0f)^2}] km/s)\n",
           Omega_p, vx_str, y0_str);
    printf("#   ND= %d  SH= %d\n", ND, SH);

    // ------- PA calculation -------
    double PA, cosPA, sinPA;
    void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
    calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);

    // ------- Sampling options -------
    auto &sampling_options = active_state->sampling;
    Initializer().read_sampling_options(context, argc, argv, cosPA, sinPA);

    // ------- Observation config -------
    ObservationConfig obs = ObservationConfig::from_cli(
        argc, argv, sampling_options.uniform_likelihood);

    // ------- Extinction -------
    if (DMrc == 0)
        DMrc = 14.3955 - 0.0239*lSIMU + 0.0122*fabs(bSIMU) + 0.128;
    if (sampling_options.n_simu == 0) exit(1);

    lDs    = (double*)malloc(sizeof(double) * 1);
    bDs    = (double*)malloc(sizeof(double) * 1);
    lDs[0] = lSIMU;
    bDs[0] = bSIMU;

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
        store_IMF_nBs(0, logMass_B, PlogM_B, PlogM_cum_norm_B, imptiles_B,
                      M0_B, M1_B, M2_B, M3_B, Ml, Mu,
                      alpha1_B, alpha2_B, alpha3_B, alpha4_B, alpha0_B);
    }

    // ------- Density grid -------
    int Dmax = getOptiond(argc, argv, "Dmax", 1, 16000);
    int npri = (ND > 0) ? 40 : 10;
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
        void calc_opticaldepth(double *tauall, double *Nsall, int idata,
                               int Dsmax21, double AI0, double hscale,
                               double Isst, double Isen);
        calc_opticaldepth(&tauall, &Nsall, 0, Dmax, AI0, hscale, Isst, Isen);
    }

    // ------- CheckD mode -------
    if (CheckD == 1) {
        double getcumu2xist(int n, double *x, double *F, double *f,
                            double Freq, int ist, int inv);
        for (int i = 0; i < 500000; i++) {
            double ran = ran1(), cumu = 0;
            int j_L, j_S;
            for (j_L = 0; j_L < ncomp; j_L++) {
                cumu += grid.cumu_rho_L_data(j_L)[grid.nbin()] /
                        grid.cumu_rho_all_L_data()[grid.nbin()];
                if (ran < cumu) break;
            }
            cumu = 0;
            for (j_S = 0; j_S < ncomp; j_S++) {
                cumu += grid.cumu_rho_S_data(j_S)[grid.nbin()] /
                        grid.cumu_rho_all_S_data()[grid.nbin()];
                if (ran < cumu) break;
            }
            double d_L = grid.sample_lens_distance(j_L, 0, grid.nbin(), ran1());
            double d_S = grid.sample_source_distance(j_S, ran1());
            void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD);
            double vxyz_L[3] = {};
            double tau_l = (j_L == 9) ? mageND + sageND*gasdev()
                         : (j_L == 8) ? mageB + sageB*gasdev()
                         : (j_L == 10) ? 14.0 : medtauds[j_L];
            get_vxyz_ran(vxyz_L, j_L, tau_l, d_L, lDs[0], bDs[0]);
            double v_l = sqrt(vxyz_L[0]*vxyz_L[0] + vxyz_L[1]*vxyz_L[1]
                            + vxyz_L[2]*vxyz_L[2]);
            printf("%d %6.0f %d %6.0f %.6f %7.2f %7.2f %7.2f %7.2f\n",
                   j_S, d_S, j_L, d_L, sqrt(d_L/8000.0),
                   vxyz_L[0], vxyz_L[1], vxyz_L[2], v_l);
        }
        exit(1);
    }

    // Release luminosity tables (no longer needed after grid build)
    population.release_luminosity_functions();

    // ------- MassFunction wrapper -------
    MassFunction mf;
    mf.init_from_population(population, context);

    // ------- Monte Carlo simulation -------
    printf("#----- Output of Monte Carlo simulation w/ VERBOSITY= %d and seed= %ld ----\n",
           sampling_options.verbosity, seed);

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
    nsd_moments.release_if_enabled();
    population.release_all();
    free(lDs);
    free(bDs);
    kinematic_tables.release_all();
    return 0;
}

} // namespace genulens
