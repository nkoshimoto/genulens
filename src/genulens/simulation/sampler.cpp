/* Generate microlensing events following the Galactic model developed by Koshimoto, Baba & Bennett (2021).
 * N. Koshimoto wrote the original .c version and C. Ranc converted it into .cpp to replace functions (ran1 and gasdev) from the Numerical Recipes in C (NR) with public alternatives (from GSL).
 * We found that the version with ran1 and gasdev from the NR was faster (~1.3 times) than the current version.
 * Please replace them by yourself if you want because we are not allowed to include a NR function in a public package. 
 * Note that a negative seed value has to be used in the NR function in contrast to a positive seed value required for this version. 
 *
 * This is version 1.2.1 of genulens.
 *
 * Version 1.2 update (Jun 14 - July 11 2022 by N.Koshimoto)
 *   1. NSD (nuclear stellar disk) component is added based on Sormani+22.
 *   2. GC is now located at the position of SgrA* and CenSgrA option is added to locate GC at (l,b) = (0,0)
 *   3. tE mean and median are calculated during the simulation (a large Nsimu value is needed to have a precision).
 *   4. CALCTAU option is added to calculate the optical depth and event rate. (event rate uses the tE mean value based on the simulation.)
 *   5. calc_PA function is added and PA is now automatically calculated from (l, b)
 *   6. Importance sampling used when tE, thetaE, and/or piE is given, which makes most calculation ~10 times faster.
 *     
 * Correct vkick bag on 2022/10/22
 * 
 * Version 1.2.1 (2023/2/18 - )
 *   Add stellar halo based on Robin+03
 * 
 * 2025/3/20
 *   Add inflation of scale height by BH kick 
 *   The same function was first installed for the analysis by Koshimoto, Kawanaka and Tsuna (2024), ApJ, 973, 5.
 *   The inflation is only for disk component, and not for bulge or other components. 
 * */
#include <math.h> 
#include <stdio.h> 
#include <unistd.h>
#include <string.h> 
#include <stdarg.h>
#include <random>
#include <memory>
#include <vector>
#include "genulens/cli/option.h"
#include <stdlib.h>
#include "genulens/io/input_data.hpp"
#include "genulens/math/interpolation.hpp"
#include "genulens/math/quadrature.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/extinction.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/rng.hpp"
#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/sampler.hpp"
#include "genulens/simulation/observation_likelihood.hpp"

#define fopen(path, mode) genulens::open_input_file((path), (mode))

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"

namespace {

struct MonteCarloStats {
  double ncntall = 0;
  double ncnts = 0;
  double ncntbWD = 0;
  double ncntbCD = 0;
  double nBD = 0;
  double nMS = 0;
  double nWD = 0;
  double nNS = 0;
  double nBH = 0;
  int Nlike = 0;
  long NrejIS = 0;
  long Ngen = 0;
  double wtlike = 0;
  double wtlike_tE = 0;
  double wtlike_except_piE = 0;
  double wtlike_w_piEe = 0;
  double wtlike_except_thE = 0;
  double wtlike_w_thEe = 0;
  double SumGamma = 0;
  double SumtE = 0;
  double logtEmin = -1;
  double logtEmax = 2;
  int NbintE = 300;
  double NlogtEs[500] = {};
};

} // namespace

#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define       PI 3.1415926535897932385
#define NDATAMAX 8000000000 // take ~6hours?
#define    KAPPA 8.1439 // 
#define STR2MIN2 8.461595e-08  // min^2 in str = deg^2 in str / 3600
#define   PI4GC2 6.013565416421e-13  // 4piG/c^2 * 1 Msun/pc = pi * kappa * mas * as
#define    KS2MY 210.949526569698696 // ([sec/yr]/[km/AU]) for km/sec/pc -> mas/yr
#define       GC 4.30091e-03 // Gravitational Constant in pc * Msun^-1 * (km/sec)^2 (Eng. Wikipedia)
#define     zsun 25.0
#define     srob 500.0
#define    vescd 550.0 // escape velo of disk
#define    vescb 600.0 // escape velo of bulge
#define   MNSMIN  1.2 // roughly determined from Fig. 7 of Raithel+18
#define   MNSMAX  2.1 // roughly determined from Fig. 7 of Raithel+18
#define  MAXMULT 1.0 // 
#define MAXGAMMA 4.0 // 
#define MINGAMMA 0.0 // 
#define MAXSIGLOGA 1.8 // 
#define MINSIGLOGA 0.3 // 
#define MAXMEANLOGA 1.7 // 
#define MINMEANLOGA 0.6 // 

namespace genulens {

genulens::RunContext *active_state = nullptr;

double ran1(){
    return active_state->runtime.rng->uniform();
}

double gasdev(){
    return active_state->runtime.rng->gaussian();
}

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

// Declare functions
double getx2y_khi(int n, double *x, double *y, double xin, int *khi);
double getx2y_ist(int n, double *x, double *y, double xin, int *ist);
double interp_x(int n, double *F, double xst, double dx, double xreq);
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq);
void   interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq);
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);

int Sampler::run_cli(RunContext &context, int argc,char **argv)
{
  active_state = &context;
  //--- Get input directory ---
  char curdir[200];
  getcwd(curdir, 200);

  //--- read parameters ---
  int CheckD   = getOptiond(argc,argv,"CheckD", 1, 0);
  Initializer().read_model_options(context, argc, argv);
  double M0_B = context.imf_options.m0;
  double M1_B = context.imf_options.m1;
  double M2_B = context.imf_options.m2;
  double M3_B = context.imf_options.m3;
  double Ml = context.imf_options.ml;
  double Mu = context.imf_options.mu;
  double alpha0_B = context.imf_options.alpha0;
  double alpha1_B = context.imf_options.alpha1;
  double alpha2_B = context.imf_options.alpha2;
  double alpha3_B = context.imf_options.alpha3;
  double alpha4_B = context.imf_options.alpha4;
  Initializer().finalize_spatial_model(context, argc, argv);

  // Read sightline and source-selection options before population setup.
  double lSIMU   = getOptiond(argc,argv,"l",  1,  1.0);
  double bSIMU   = getOptiond(argc,argv,"b",  1, -3.9);
  double Isst    = getOptiond(argc,argv,"Isrange", 1, 14.0);
  double Isen    = getOptiond(argc,argv,"Isrange", 2, 21.0);
  double VIsst   = getOptiond(argc,argv,"VIsrange", 1, 0.0);
  double VIsen   = getOptiond(argc,argv,"VIsrange", 2, 0.0);
  double AIrc    = getOptiond(argc,argv,"AIrc", 1, 0);
  double EVIrc   = getOptiond(argc,argv,"EVIrc", 1, 0);
  double DMrc    = getOptiond(argc,argv,"DMrc", 1, 0);
  double AKrc    = getOptiond(argc,argv,"AKrc", 1, 0);
  if (fabs(lSIMU) > 10 || fabs(bSIMU) > 7 || fabs(bSIMU) < 1.5){
    printf("# WARNING: genulens is designed for use in |l| < ~10 and ~1.5 < |b| < ~7 and (l,b)= (%.3f, %.3f) deg. is outside of the range.\n",lSIMU,bSIMU);
  }

  PopulationRuntime population;
  population.initialize_mass_function(context.imf_options);
  population.initialize_luminosity_functions(Isst, Isen, VIsst, VIsen, AIrc, EVIrc);
  population.read_empirical_mass_luminosity();
  double *logMass_B = population.log_mass;
  double *PlogM_B = population.mass_probability;
  double *PlogM_cum_norm_B = population.mass_cumulative;
  int *imptiles_B = population.mass_percentiles;
  double *M_emps = population.empirical_masses;
  double **Mag_emps = population.empirical_magnitudes;
  int nMLemp = population.empirical_count;

  KinematicRuntimeTables kinematic_tables;
  kinematic_tables.initialize_shu_distribution();

  // set y0d for disk normalize
  y0d[0] = (DISK == 1) ? exp(-R0/Rd[0] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[0]);
  y0d[1] = (DISK == 1) ? exp(-R0/Rd[1] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[1]);
  y0d[2] = (DISK == 1) ? exp(-R0/Rd[2] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[2]);

  /*------------ For BH kick analysis --------------------*/
  //---- Set kick velocity  -----
  int    MXDkick = getOptiond(argc,argv,"MXDkick",  1,    0); //  Give v_kick in maxwell distribution
  double vkickBH = getOptiond(argc,argv,"vkickBH",  1,  100); //  Table 2 of Lam et al. 2020
  double vkickNS = getOptiond(argc,argv,"vkickNS",  1,  350); //  Table 2 of Lam et al. 2020
 
  //------ Parameters to change scale height for BH (only for disk as of 2022/10/21)
  //  hd(R,vkick) = hd0 * exp[R/Rhd] * (vkick/[50km/s])^beta_v
  //  hd(R,vkick) = hd0 * exp[R/Rhd] * (vavg/[50km/s])^beta_v <- should be used?
  //  Values from ../BHkick/anal_kick_velo_20221027.ipynb based on Tsuna+2018
  int    BHhd  = getOptiond(argc,argv,"BHhd",  1,    0); // flag to change BH scale height for disk
  int    BHhb  = getOptiond(argc,argv,"BHhb",  1,    0); // flag to change BH scale height for bar (just naively)
  int fixRhdBH = getOptioni(argc,argv,"fixRhdBH",  1,    0); // flag to change BH scale height for bar (just naively)
  double RhdBH0 = getOptiond(argc,argv,"RhdBH0",  1, 9660); // scale length [pc] of hd(R)
  double betaBH = getOptiond(argc,argv,"betaBH",  1, 0.820); // power of vkick dependence.
  /*------------------------------------------------------*/

  //------ Parameters to change surface density (added on 2022/11/23, see anal_kick_velo_20221123.ipynb)
  //  Sigma(R, vkick) = Sigma(R, 50km/s) * (a2*R^2 + a1*R + a0)
  //  a2, a1, a0 are determined for vkick = 100, 200, 400 km/s in anal_kick_velo_20221123.ipynb
  int  UseSigBH = getOptiond(argc,argv,"UseSigBH",  1,    0); // Use Sigma dependency on vkick
  double a2toSig25BHs[5] = {0, 0.00279,  0.00985,  0.02054,  0.01901}; // coeff of R^2 of Sig_v / Sig_25
  double a1toSig25BHs[5] = {0,-0.02023, -0.07548, -0.16585, -0.15082}; // coeff of R^1 of Sig_v / Sig_25
  double a0toSig25BHs[5] = {1, 1.01342,  1.05382,  1.05898,  0.69889}; // coeff of R^0 of Sig_v / Sig_25

  int  printBHfac = getOptiond(argc,argv,"printBHfac",  1,    0); // print fBH 

  // Print input parameters as header 
  printf("#   Output of \"./genulens ");
  for (int i=1;i<argc;i++) {
    printf("%s", argv[i]);
    if (i < argc - 1) {
      printf(" ");
    }
  }
  printf("\"\n");
  printf("#   You need to weight each line by the 1st value of each line from the left \n");
  printf("#---------- Parameters for IMF and Sun ----------\n");
  printf("#      IMF:  alpha0= %5.2f ( %.2f <M< %.2f ),\n",alpha0_B,M0_B,Mu);
  printf("#            alpha1= %5.2f ( %.2f <M< %.2f ),\n",alpha1_B,M1_B,M0_B); 
  printf("#            alpha2= %5.2f ( %.2f <M< %.2f ),\n",alpha2_B,M2_B,M1_B); 
  printf("#            alpha3= %5.2f ( %.2f <M< %.2f ),\n",alpha3_B,M3_B,M2_B); 
  printf("#            alpha4= %5.2f ( %.5f <M< %.5f )\n",alpha4_B,Ml,M3_B);
  printf("#         (R, z)sun= (%5.0f, %5.0f) pc\n",R0,zsun);
  printf("#   (vx, vy, vz)sun= (%5.1f, %5.1f, %4.1f) km/s\n",vxsun,vysun,vzsun);
  printf("#------------ Disk model: (DISK, hDISK, tSFR)= ( %d , %d , %.1f Gyr ) --------------\n",DISK, hDISK,tSFR);
  printf("#            tau   Rd  zd zd45 sigU0 sigW0  RsigU  RsigW    rho0        n0     n0WD    Sigma0  \n");
  printf("#            Gyr   pc  pc   pc  km/s  km/s     pc     pc  Msun/pc^3  */pc^3   */pc^3  Msun/pc^2\n");
  double MVVVd = 0;  // mass in the VVV box when DISK == 2
  double Mind = 0;  // mass of inner disk (< Rbreak) when DISK == 2
  for (int i = 0; i< 8; i++){
    int rd = (i == 0) ? Rd[0] : (i < 7) ? Rd[1] : (i == 7) ? Rd[2] : 0;
    int zdtmp = (hDISK == 0) ? zd[i] : zd45[i];
    double MVVVtmp = 0;
    // if (DISK == 2){ // same normalization also when DISK != 2
      MVVVtmp = rho0d[i] * exp((R0 - Rdbreak)/rd) * 2200*2 * 1400*2 * zd[i] / zdtmp;
      double ztmp    = 1200.0/zdtmp;  // z of VVV box
      Mind    += 2 * zdtmp * MVVVtmp / 4400 / 2800 * PI * Rdbreak * Rdbreak;
      MVVVtmp *= (i < 7) ? 2 * zdtmp * (exp(2*ztmp) - 1)/(exp(2*ztmp) + 1)  
                         : 2 * zdtmp * (1 - exp(-ztmp));
      MVVVd   += MVVVtmp;
    // }
    double hsigU = (i < 7) ? hsigUt : hsigUT;
    double hsigW = (i < 7) ? hsigWt : hsigWT;
    double sigW0 = (i < 7) ? sigW10d * pow((medtauds[i]+0.01)/10.01, betaW) : sigW0td;
    double sigU0 = (i < 7) ? sigU10d * pow((medtauds[i]+0.01)/10.01, betaU) : sigU0td;
    double SigmaD = 2 * zd[i] * rho0d[i];
    printf ("#   Disk%d: %5.2f %4d %3.0f  %3d %5.2f %5.2f %6.0f %6.0f   %.2e %.2e %.2e   %.2e\n",i+1,medtauds[i], rd, zd[i],zdtmp,sigU0,sigW0,hsigU,hsigW,rho0d[i],n0d[i],n0d[i]-n0MSd[i], SigmaD);
  }
  printf("#------------ Disk BH scale height --------------\n");
  printf("#            BHhd      = %d (0: hdBH = hdstar, 1: hd(R,vkick) = hd0*exp[R/Rhd]*(vkick/[50km/s])^beta_v)\n", BHhd);
  if (BHhd == 1){
    printf ("#   (Rhd, beta_v) =  (%6.1f pc, %6.4f) \n", RhdBH0, betaBH);
  }
  printf ("#        UseSigBH      = %d (0: Sigma(R) same against vkick, 1: Sigma(R) depends on vkick\n", UseSigBH);
  printf ("#         MXDkick      = %d (0: vkick constant, 1: vkick given by Maxwell distribution)\n", MXDkick);
  printf ("#   (vkickNS, vkickBH) =  (%6.1f , %6.1f) km/s\n", vkickNS, vkickBH);

  // Crude normalize bulge mass
  double crude_integrate(double xmax, double ymax, double zmax, int nbun);
  double massVVVbox = crude_integrate(2200, 1400, 1200, 15); // VVV box defined by Wegg & Gerhard (2013), MNRAS, 435, 1874
  double massentire = crude_integrate(6000, 3000, 3000, 30); // should include entire bulge
  double fm1 = 1, fmX = 0;
  if (addX >= 5){
    int addXtmp = addX;
    addX = 0;
    double mass1all = crude_integrate(6000, 3000, 3000, 30);
    addX = addXtmp;
    fm1 = mass1all / massentire;
    fmX = 1 - fm1;
  }  
  double MVVVP17 = 1.32e+10;
  rho0b = (frho0b * MVVVP17 - MVVVd)/massVVVbox; // normalized by mass in the VVV box by Portail et al. (2017), MNRAS, 465, 1621
  n0MSb = rho0b * fb_MS * m2nb_MS; // number density of bar MS stars
  n0RGb = n0MSb * nMS2nRGb; // number density of bar RG stars (for mu calculation)
  n0b   = n0MSb + rho0b * (1 - fb_MS) * m2nb_WD; // number density of b MS+WD stars
  massVVVbox *= rho0b;
  massentire *= rho0b;

  // Consider Nuclear Disk if  (y, z) reaches (125, 50) x 5 (= 625, 250) at 8 kpc
  double MND;
  if (fabs(lSIMU) < 5 && fabs(bSIMU) < 2)  ND = 3;
  ND       = getOptiond(argc,argv,"NSD",     1,  ND); // 0: wo nuclear disk, 1: w/ nuclear disk by P17, 2: w/ Sormani+21-like NSD
  if (ND > 0) ND = 3;
  if (ND == 1){ // Consider Portail+17's NSD
    MND  = 2.0e+09;
    x0ND = 250;
    y0ND = 125;
    z0ND =  50;
  }
  if (ND == 2){ // Consider Sormani+21-like NSD
    MND  = 7.0e+08;
    x0ND =  74;
    y0ND =  74;
    z0ND =  26;
  }
  x0ND  = getOptiond(argc,argv,"x0ND",  1, x0ND); 
  y0ND  = getOptiond(argc,argv,"y0ND",  1, y0ND); 
  z0ND  = getOptiond(argc,argv,"z0ND",  1, z0ND); 
  MND   = getOptiond(argc,argv,"MND" ,  1,  MND); 
  // normalize ND mass
  if (ND){
    rho0ND = (ND == 3) ? 1 : 0.25*MND/PI/x0ND/y0ND/z0ND; // Msun/pc^3 is given by calc_rho_each when ND == 3
    n0MSND = rho0ND * fND_MS * m2nND_MS; // number density of ND MS stars
    n0RGND = n0MSND * nMS2nRGND; // number density of ND RG stars (for mu calculation)
    n0ND   = n0MSND + rho0ND * (1 - fND_MS) * m2nND_WD; // number density of ND MS+WD stars
  }
  nzND = (zenND - zstND)/dzND + 1.5;
  nRND = (RenND - RstND)/dRND + 1.5;
  if (ND == 3){ // More Sormani+21-like NSD, Use input_files/NSD_moments.dat 
    logrhoNDs   = (double**)malloc(sizeof(double *) * nzND);
    vphiNDs     = (double**)malloc(sizeof(double *) * nzND);
    corRzNDs    = (double**)malloc(sizeof(double *) * nzND);
    logsigvNDs  = (double***)malloc(sizeof(double *) * nzND);
    for (int i=0; i<nzND; i++){
      logrhoNDs[i] = (double*)calloc(nRND, sizeof(double *));
      vphiNDs[i]   = (double*)calloc(nRND, sizeof(double *));
      corRzNDs[i]  = (double*)calloc(nRND, sizeof(double *));
      logsigvNDs[i] = (double**)malloc(sizeof(double *) * nRND);
      for (int j=0; j<nRND; j++){
        logsigvNDs[i][j] = (double*)calloc(3, sizeof(double *)); // 3= phi, R, z
      }
    }
    char *fileND = (char*)"input_files/NSD_moments.dat";
    void store_NSDmoments(char *infile);
    store_NSDmoments(fileND);
  }
  printf("#------------------ Bulge model: (alpha_bar, Mbar, Mind, MVVVb, MVVVd) = ( %.1f deg, %.2e Msun, %.2e Msun, %.2e Msun, %.2e Msun) ---------------------\n",thetaD,massentire,Mind,massVVVbox,MVVVd);
  printf("#   (M_MS, M_REM)ave= (%.6f %.6f) Msun/*, fM_REM= %.4f, Mass/RG= %5.1f Msun/RG \n",1/m2nb_MS,1/m2nb_WD,1-fb_MS,1/fb_MS/m2nb_MS/nMS2nRGb);
  printf("#   rho%d: M= %.2e Msun, rho0b= %5.2f Msun/pc^3, (x0, y0, z0, Rc)= (%4.0f, %4.0f, %3.0f, %4.0f) pc, (C1, C2,   C3)= (%.1f, %.1f, %.1f)\n",model,fm1*massentire,rho0b,x0_1,y0_1,z0_1,Rc,C1,C2,C3);
  if (addX >= 5) printf("#     X%d: M= %.2e Msun, rho0X= %5.2f Msun/pc^3, (x0, y0, z0, Rc)= (%4.0f, %4.0f, %3.0f, %4.0f) pc, (C1, C2, b_zX)= (%.1f, %.1f, %.1f)\n",addX,fmX*massentire,rho0b*fX,x0_X,y0_X,z0_X,Rc_X,C1_X,C2_X,b_zX);
  printf("#   (Omega_p, vx_str)= ( %.1f km/s/kpc, %3.0f[1 - e^{-(|yb|/%4.0f)^2}] km/s ),",Omega_p,vx_str,y0_str);
  printf(" sig0+1(xb, yb, zb)= (%3.0f+%3.0f, %3.0f+%3.0f, %3.0f+%3.0f) km/s\n",sigx_vb,sigx_vb0,sigy_vb,sigy_vb0,sigz_vb,sigz_vb0);
  printf("#   sigR%d: (x0, y0, z0)= (%5.0f, %5.0f, %5.0f) pc, (C1, C2, C3)= (%.1f, %.1f, %.1f)\n",model_vb,x0_vb,y0_vb,z0_vb,C1_vb,C2_vb,C3_vb);
  printf("#   sigZ%d: (x0, y0, z0)= (%5.0f, %5.0f, %5.0f) pc, (C1, C2, C3)= (%.1f, %.1f, %.1f)\n",model_vbz,x0_vbz,y0_vbz,z0_vbz,C1_vbz,C2_vbz,C3_vbz);
  printf("#   ND= %d     (0: no NSD, 1: Portail+17's NSD, 2: Sormani+22-like NSD, 3: Use Sormani+22's DF's moments)\n", ND);
  printf("#   SH= %d     (0: no stellar halo, 1: Robin+03's Spheroid)\n", SH);
  // printf("#   n0MSSH= %.5e, n0SH= %.5e\n",n0MSSH, n0SH);
  
  // Calc PA (added on 2022/6/14)
  double PA, cosPA, sinPA;
  void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
  calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);
  
  // Read Input parameters for simulation
  auto &sampling_options = active_state->sampling;
  Initializer().read_sampling_options(context, argc, argv, cosPA, sinPA);
#define NSIMU sampling_options.n_simu
#define NlikeMIN sampling_options.n_like_min
  // double PA         = getOptiond(argc,argv,"PA", 1, 59.56); // postition angle l to E value. 59.56 is for (l,b) = (0.94, -1.48), needed to calculate piEn and piEe
  // double cosPA    = cos(PA/180.0*PI), sinPA = sin(PA/180.0*PI);
  // printf("PA= %f, cosPA= %f, sinPA= %f\n",PA,cosPA,sinPA);
#define vEarthl sampling_options.v_earth_l
#define vEarthb sampling_options.v_earth_b
#define gammaDs sampling_options.gamma_ds
#define wtD_L sampling_options.weight_lens_distance
#define wtM_L sampling_options.weight_lens_mass
#define NoGAMMAIS sampling_options.no_gamma_importance_sampling
#define SMALLGAMMA sampling_options.small_gamma
#define VERBOSITY sampling_options.verbosity
#define UNIFORM sampling_options.uniform_likelihood
#define BINARY sampling_options.binary
#define REMNANT sampling_options.remnant
#define onlyWD sampling_options.only_white_dwarf
#define CALCPRIORpiE sampling_options.calc_prior_piE
#define CALCPRIORthE sampling_options.calc_prior_thetaE
  double tEobs     = getOptiond(argc,argv,"tE", 1, 54.5); // in day
  double tEe       = getOptiond(argc,argv,"tE", 2, 99999999999.0); // in day
  double fetE = 0;
  int    tEdet     = getOptiond(argc,argv,"tEdet", 1, 0); // 0: det, 1: upper limit, 2: lower limit
  double tEmin = getOptiond(argc,argv,"tErange", 1, 0); // in mas
  double tEmax = getOptiond(argc,argv,"tErange", 2, 0); // in mas
  double thetaEobs = getOptiond(argc,argv,"thetaE", 1, 0); // in mas
  double thetaEe   = getOptiond(argc,argv,"thetaE", 2, 0); // in mas
  double fethetaE = 0;
  int    thetaEdet = getOptiond(argc,argv,"thetaEdet", 1, 0); // 0: det, 1: upper limit, 2: lower limit
  double thetaEmin = getOptiond(argc,argv,"thetaErange", 1, 0); // in mas
  double thetaEmax = getOptiond(argc,argv,"thetaErange", 2, 0); // in mas
  double piEobs = getOptiond(argc,argv,"piE", 1, 0); // parallax amplitude 
  double piEe   = getOptiond(argc,argv,"piE", 2, 0); // 
  double fepiE = 0;
  int    piEdet = getOptiond(argc,argv,"piEdet", 1, 0); // 0: det, 1: upper limit, 2: lower limit
  double piEmin = getOptiond(argc,argv,"piErange", 1, 0); 
  double piEmax = getOptiond(argc,argv,"piErange", 2, 0); 
  double piENobs = getOptiond(argc,argv,"piEN", 1, 0); // 
  double piENe   = getOptiond(argc,argv,"piEN", 2, 0); // 
  double fepiEN = 0;
  double piEEobs = getOptiond(argc,argv,"piEE", 1, 0); // 
  double piEEe   = getOptiond(argc,argv,"piEE", 2, 0); // 
  double fepiEE = 0;
  double muslobs = getOptiond(argc,argv,"musl", 1, 0); // in mas/yr 
  double musle   = getOptiond(argc,argv,"musl", 2, 0); // in mas/yr 
  double femusl = 0;
  double musbobs = getOptiond(argc,argv,"musb", 1, 0); // in mas/yr
  double musbe   = getOptiond(argc,argv,"musb", 2, 0); // in mas/yr
  double femusb = 0;
  double musNobs = getOptiond(argc,argv,"musN", 1, 0); // in mas/yr
  double musNe   = getOptiond(argc,argv,"musN", 2, 0); // in mas/yr
  double femusN = 0;
  double musEobs = getOptiond(argc,argv,"musE", 1, 0); // in mas/yr
  double musEe   = getOptiond(argc,argv,"musE", 2, 0); // in mas/yr
  double femusE = 0;
  double muhelNobs = getOptiond(argc,argv,"muhelN", 1, 0); // in mas/yr
  double muhelNe   = getOptiond(argc,argv,"muhelN", 2, 0); // in mas/yr
  double femuhelN  = 0; // 
  double muhelEobs = getOptiond(argc,argv,"muhelE", 1, 0); // in mas/yr
  double muhelEe   = getOptiond(argc,argv,"muhelE", 2, 0); // in mas/yr
  double femuhelE  = 0; // 
  int    musRCG  = getOptiond(argc,argv,"musRCG", 1, 0); // 0: mus relative to Sun, 1: mus relative to RCG 
  double ILobs = getOptiond(argc,argv,"IL", 1, 14.00); // mag
  double ILe   = getOptiond(argc,argv,"IL", 2,  0.01); // mag
  // double feIL  = getOptiond(argc,argv,"IL", 3, 0); //
  double feIL = 0; 
  int    ILdet = getOptiond(argc,argv,"ILdet", 1, 2); // 0: det, 1: upper limit, 2: lower limit, default: 2
  double KLobs = getOptiond(argc,argv,"KL", 1,  0); // mag
  double KLe   = getOptiond(argc,argv,"KL", 2,  0); // mag
  // double feKL  = getOptiond(argc,argv,"KL", 3, 0); //
  double feKL = 0; 
  int    KLdet = getOptiond(argc,argv,"KLdet", 1, 0); // 0: det, 1: upper limit, 2: lower limit, default: 2
  double u0obs = getOptiond(argc,argv,"u0", 1, 0); // affect only when BINARY == 1
  // char *pthEfile = getOptions(argc,argv, "pthetaE", 1, ""); // give P(thE) by file
  //
  // Set importance sampling parameters
  int    NOIS = getOptiond(argc,argv,"NOIS", 1, 0); // 0: Set tErange, thetaErange, piErange automatically
  double fIS0 = (UNIFORM == 1) ? 1.02 : 4.0;
  double fIStE     = getOptiond(argc,argv,"fIStE", 1, fIS0); // Consider +- 4 sigma
  double fISthetaE = getOptiond(argc,argv,"fISthetaE", 1, fIS0); // Consider +- 4 sigma
  double fISpiE    = getOptiond(argc,argv,"fISpiE", 1, fIS0); //
  if (NOIS == 0){
    if (tEobs - tEe > 0 && tEe > 0 && tEmax - tEmin == 0){
      tEmin = tEobs - fIStE * tEe;
      tEmax = tEobs + fIStE * tEe;
      if (tEmin <= 0 || tEdet == 1) tEmin = 1e-10;
      if (tEdet == 2) tEmax = 1e+6;
    }
    if (thetaEobs - thetaEe > 0 && thetaEe > 0 && thetaEmax - thetaEmin == 0){
      thetaEmin = thetaEobs - fISthetaE * thetaEe;
      thetaEmax = thetaEobs + fISthetaE * thetaEe;
      if (thetaEmin <= 0 || thetaEdet == 1) thetaEmin = 1e-10;
      if (thetaEdet == 2) thetaEmax = 1e+6;
    }
    if (piEobs - piEe > 0 && piEe > 0 && piEmax - piEmin == 0){
      piEmin = piEobs - fISpiE * piEe;
      piEmax = piEobs + fISpiE * piEe;
      if (piEmin <= 0 || piEdet == 1) piEmin = 1e-10;
      if (piEdet == 2) piEmax = 1e+6;
    }
    if (piENe > 0 && piEEe > 0 && piEmax - piEmin == 0){
      double piENabsmin = (fabs(piENobs) - fISpiE * piENe < 0) ? 0 : fabs(piENobs) - fISpiE * piENe;
      double piEEabsmin = (fabs(piEEobs) - fISpiE * piEEe < 0) ? 0 : fabs(piEEobs) - fISpiE * piEEe;
      double piENabsmax = fabs(piENobs) + fISpiE * piENe;
      double piEEabsmax = fabs(piEEobs) + fISpiE * piEEe;
      piEmin = sqrt(piENabsmin*piENabsmin + piEEabsmin*piEEabsmin);
      piEmax = sqrt(piENabsmax*piENabsmax + piEEabsmax*piEEabsmax);
      if (piEmin <= 0) piEmin = 1e-10;
    }
  }

  if (DMrc == 0)
    DMrc = 14.3955 - 0.0239 * lSIMU + 0.0122*fabs(bSIMU)+0.128; // Eqs(2)-(3) of Nataf+16 

  if (NSIMU == 0) exit(1);
  int idata = 0;
  lDs        = (double *)malloc(sizeof(double *) * 1);
  bDs        = (double *)malloc(sizeof(double *) * 1);
  lDs[0]    = lSIMU; //
  bDs[0]    = bSIMU; //
  double hdust = getOptiond(argc,argv,"hdust", 1, 164.0); // in pc. 164 pc = dust scale height from Nataf+13 
  double cosb = cos(bDs[idata]/180.0*PI), sinb = sin(bDs[idata]/180.0*PI), 
         cosl = cos(lDs[idata]/180.0*PI), sinl = sin(lDs[idata]/180.0*PI);
  double hscale = hdust/(fabs(sinb) + 0.0001);  // 164 pc = dust scale height from Nataf+13
  double Dmean  = (DMrc > 0) ? pow(10, 0.2*DMrc) * 10
                 : -9.99;
  gmodel::ExponentialDustExtinction extinction(hscale, Dmean, AIrc, AKrc, EVIrc);
  double AI0 = extinction.ai0();
  double AK0 = extinction.ak0();
  double EVI0 = extinction.evi0();
  // if (AI0 == 0){ // not recommended, AIrc should be given if known
  //   AI0 = 0.64 * hscales[0] * 0.001 -0.33; // by linear fit to known AIs of Nataf+13
  //   if (AI0 > 7) AI0 = 7;
  //   if (AI0 < 0) AI0 = 0;
  // }
  // if ((thetaEmax - thetaEmin > 0 || piEmax - piEmin > 0) && SMALLGAMMA == 0){
  //   printf("# SMALLGAMMA set to 1 because thetaErange or piErange is given\n");
  //   SMALLGAMMA = 1;
  // }
  printf("#-------------- Input parameters ---------------\n");
  printf("#    CenSgrA= %d     (0: GC at (l,b)=(0,0), 1: GC at (l,b)= (%.3f, %.3f))\n",
         context.spatial.center_on_sgr_a, context.spatial.l_sgr_a, context.spatial.b_sgr_a);
  printf("#    UNIFORM= %d     (0: L= N(obs, err), 1: L= U(obs-err, obs+err))\n", UNIFORM);
  printf("#    REMNANT= %d     (0: no remnant, 1: with remnant)\n", REMNANT);
  printf("#     onlyWD= %d     (1: w/ remnant but only WD)\n", onlyWD);
  printf("#     BINARY= %d     (0: no binary , 1: with binary )\n", BINARY);
  printf("# SMALLGAMMA= %d     (1: output all events even with small Gamma )\n", SMALLGAMMA);
  printf("#   (Nsimu, NlikeMIN)= (%ld, %ld)\n",NSIMU,NlikeMIN);
  printf("#          (l, b, PA)= (%6.2f, %6.2f, %5.2f) deg.\n",lSIMU,bSIMU,PA);
  printf("#       (vl, vb)Earth= (%6.2f, %6.2f) km/s\n",vEarthl,vEarthb);

  if (tEe > 0)     printf("#     tE = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",tEobs,tEe,fetE,tEdet);
  if (thetaEe > 0) printf("# thetaE = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",thetaEobs,thetaEe,fethetaE,thetaEdet);
  if (piEe > 0)    printf("#    piE = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",piEobs,piEe,fepiE,piEdet);
  if (piENe > 0 && piEEe > 0) 
     printf("#    (piEN, piEE) = (%.3f, %.3f) +- (%.3f, %.3f) Fe_add= (%6.1f, %6.1f)\n",piENobs,piEEobs,piENe,piEEe,fepiEN,fepiEE);
  if (musle > 0 && musbe > 0) 
     printf("#    (musl, musb) = (%.3f, %.3f) +- (%.3f, %.3f) Fe_add= (%6.1f, %6.1f)\n",muslobs,musbobs,musle,musbe,femusl,femusb);
  if (musNe > 0 && musEe > 0) 
     printf("#    (musN, musE) = (%.3f, %.3f) +- (%.3f, %.3f) Fe_add= (%6.1f, %6.1f)\n",musNobs,musEobs,musNe,musEe,femusN,femusE);
  if (muhelNe > 0 && muhelEe > 0) 
     printf("#    (muhelN, muhelE) = (%.3f, %.3f) +- (%.3f, %.3f) Fe_add= (%6.1f, %6.1f)\n",muhelNobs,muhelEobs,muhelNe,muhelEe,femuhelN,femuhelE);
  if (ILe > 0)    printf("#      IL = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",ILobs,ILe,feIL,ILdet);
  if (KLe > 0)    printf("#      KL = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",KLobs,KLe,feKL,KLdet);
  if (AK0   > 0)   printf("#     Consider (hdust, Dmean, AKrc,  AK0) = (%.0f, %.0f, %.2f, %.2f)\n",hdust,Dmean,AKrc,AK0);
  if (u0obs > 0)   printf("#  u0obs = %.3f\n",u0obs);
  if (AI0   > 0)   printf("#  Consider %.2f <  Is < %.2f, (hdust, Dmean,  AIrc,  AI0) = (%.0f, %.0f, %.2f, %.2f)\n",Isst,Isen,hdust,Dmean,AIrc,AI0);
  if (EVI0  > 0)   printf("#  Consider %.2f < VIs < %.2f, (hdust, Dmean, EVIrc, EVI0) = (%.0f, %.0f, %.2f, %.2f)\n",VIsst,VIsen,hdust,Dmean,EVIrc,EVI0);
  if (AI0 == 0 && EVI0 == 0) printf ("# gammaDs=    %.2f      : omomi in Gamma as Ds^gammaDs\n",gammaDs);
  if (tEmax - tEmin > 0) 
    printf("# Sampling parameter range : \n");
    printf("#     tErange     : %.4f - %.4f\n",tEmin, tEmax);
  if (thetaEmax - thetaEmin > 0) 
    printf("#     thetaErange : %.4f - %.4f\n",thetaEmin, thetaEmax);
  if (piEmax - piEmin > 0) 
    printf("#     piErange    : %.4f - %.4f\n",piEmin, piEmax);


  // Weight by wtM_L
  if (wtM_L != 0){
    alpha0_B += wtM_L;
    alpha1_B += wtM_L;
    alpha2_B += wtM_L;
    alpha3_B += wtM_L;
    alpha4_B += wtM_L;
    store_IMF_nBs(0, logMass_B, PlogM_B, PlogM_cum_norm_B, imptiles_B, M0_B, M1_B, M2_B, M3_B, Ml, Mu, alpha1_B, alpha2_B, alpha3_B, alpha4_B, alpha0_B);
  }



  //------- Store cumu_rho for each ith disk as a function of distance -----------
  void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb);  // return rho for each component 
  double rhos[11] = {}, xyz[3] = {}, xyb[2] = {}, sigvbs[3] = {};
  // Lens   : include REMNANT, mass basis 
  // Source : only stars, number basis 
  // int Dmax = 16000;
  int Dmax = getOptiond(argc,argv,"Dmax",  1,  16000); // Does not work due to Shu's DF
  // int Dmax = 12000;
  int nbin = (ND > 0 && fabs(lSIMU) < 0.05 && fabs(bSIMU) < 0.05) ? 0.20*Dmax+0.5 
           : (ND > 0 && fabs(lSIMU) < 0.10 && fabs(bSIMU) < 0.10) ? 0.10*Dmax+0.5
           : (ND > 0) ? 0.04*Dmax+0.5 
           : 0.01*Dmax+0.5;
  double dD = (double) Dmax/nbin;
  double *D, **rhoD_S, **rhoD_L, **cumu_rho_S, **cumu_rho_L, *cumu_rho_all_S, *cumu_rho_all_L, **fBH;
  D               = (double *)calloc(nbin+1, sizeof(double *));
  cumu_rho_all_S  = (double *)calloc(nbin+1, sizeof(double *));
  cumu_rho_all_L  = (double *)calloc(nbin+1, sizeof(double *));
  fBH         = (double **)malloc(sizeof(double *) * ncomp);  // store rho_BH/rho_*
  rhoD_S      = (double **)malloc(sizeof(double *) * ncomp);
  rhoD_L      = (double **)malloc(sizeof(double *) * ncomp);
  cumu_rho_S  = (double **)malloc(sizeof(double *) * ncomp);
  cumu_rho_L  = (double **)malloc(sizeof(double *) * ncomp);
  for (int i=0; i<ncomp; i++){
    fBH[i]        = (double *)calloc(nbin+1, sizeof(double *));
    rhoD_S[i]     = (double *)calloc(nbin+1, sizeof(double *));
    rhoD_L[i]     = (double *)calloc(nbin+1, sizeof(double *));
    cumu_rho_S[i] = (double *)calloc(nbin+1, sizeof(double *));
    cumu_rho_L[i] = (double *)calloc(nbin+1, sizeof(double *));
  }
  double fLF_detect(double extI, double Imin, double Imax, int idisk);
  double fIVI_detect(double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk);
  printf("#----- Mass density distribution along (l, b)= (%.3f, %.3f) w/ wtD_L= %.1f --------\n",lSIMU,bSIMU,wtD_L);
  int npri = (ND > 0) ? 40 : 10;
  if (printBHfac)
    npri = 100000;
  npri = getOptioni(argc,argv,"npri",  1,  npri); 
  int printrhoS = getOptiond(argc,argv,"printrhoS",  1,  0);
  double nallS = 0;
  for (int ibin=0; ibin<=nbin; ibin++){
    D[ibin] = (double) ibin/nbin * Dmax;
    calc_rho_each(D[ibin], idata, rhos, xyz, xyb);
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double z = xyz[2]; 
    if (ibin%npri==0) printf ("# %5.0f %5.0f %5.0f ",D[ibin],R,z);
    double rhosum = 0;
    const auto dust = extinction.at_distance(D[ibin]);
    double extI  =  dust.i_band + extinction.distance_modulus_term(D[ibin] + 0.1);
    double extVI = dust.color_vi;
    // printf ("%5.0f %7.3f ",D[ibin],extI);
    for (int i=0;i<ncomp;i++){
      double fBHtmp = 1;
      if (UseSigBH == 1 && i < 9 && vkickBH > 25){
        // Change BH surface density following Eq. (10) of Koshimoto+24
        if (vkickBH != 25 && vkickBH != 50 && vkickBH != 100 && vkickBH != 200 && vkickBH != 400){
          printf("# WARNING : vkickBH = 25, 50, 100, 200, or 400 is preferable for UseSigBH == 1\n");
        }
        int ivkick = (int) floor(log(vkickBH/25)/log(2) + 0.5);
        double corSigFac = a2toSig25BHs[ivkick]*R*R*1e-06 + a1toSig25BHs[ivkick]*R*1e-03 + a0toSig25BHs[ivkick];
        fBHtmp *= corSigFac;
      }
      int ien = (BHhb == 1) ? 9 : 8;
      if (BHhd == 1 && i < ien){ // only for Disk when BHhb == 0, else also for Bar (crude attempt though)
        // Change BH scale height (disk: Eq. 5 of Koshimoto+24, bar: Eq. 11 of Koshimoto+24)
        double sigW0 = (i < 7) ? sigW10d * pow((medtauds[i]+0.01)/10.01, betaW) : sigW0td;
        double hsigW = (i < 7) ? hsigWt : hsigWT;
        if (BHhb == 1 && i == 8){
           void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
           calc_sigvb(xyb[0], xyb[1], z,  sigvbs);
        }
        double sigW  = (i == 8) ? sigvbs[2] : sigW0*exp(-(R - R0)/hsigW);
        double sigW2 = sigW*sigW;
        double sigzkick2 = vkickBH*vkickBH*PI/8; // v_avg = 2 sigma_1D sqrt(2/pi)
        double sigzadd = sqrt(sigW2 + sigzkick2);
        double fvBH = pow((sigzadd/sigW), betaBH);
        double RhdBH = (fixRhdBH == 1) ? RhdBH0 : RhdBH0 * (1 + sigW2/sigzkick2); // Eq.(6) of the paper
        double fzdBH = exp((R-R0)/RhdBH) * fvBH;
        double rhoMS = 0, rhoBH = 0;
        double zd0  = (i < 8) ? zd[i] : 235.344943180979;  // use z0_1 of E model for bar.
        double zdBH = (i < 8) ? zd0 * fzdBH : zd0 * fvBH;  // doesn't consider the exp((R-R0)/RhdBH) term for bar 
        if (i == 8){ // bar, calc rhoB for E_fg0 model, when BHhb == 1
          double x0E = 668.323640191308, y0E = 277.674592258175; 
          double C1E = 1.40903573470129, C2E = 3.3497118832179;
          double xn = fabs(xyb[0]/x0E), yn = fabs(xyb[1]/y0E);
          double Rs = pow((pow(xn, C1E) + pow(yn, C1E)), 1/C1E);
          double zn, rs;
          zn = fabs(z/zd0);
          rs = pow(pow(Rs, C2E)  + pow(zn, C2E), 1/C2E);
          rhoMS = exp(-rs);
          /* Calc SigZMS here if done? */
          zn = fabs(z/zdBH);
          rs = pow(pow(Rs, C2E)  + pow(zn, C2E), 1/C2E);
          rhoBH = exp(-rs) / fvBH; // divided by fvBH to reduce rho0 (cuz Sigma is still conserved at this stage)
          /* Calc SigZBH here if done? */
        }else{
          rhoMS = (i < 7) ? 4.0/(exp(2*z/zd0)+exp(-2*z/zd0)+2)
                          : exp(-fabs(z)/zd0);
          rhoBH = (i < 7) ? 4.0/(exp(2*z/zdBH)+exp(-2*z/zdBH)+2)
                          : exp(-fabs(z)/zdBH);
          rhoBH /= fzdBH; // divided by fzdBH to reduce rho0 (cuz Sigma is still conserved at this stage)
        }
        if (printBHfac)
          printf("%d %5.0f %4.0f %5.2f %6.2f %6.2f %.6f %6.1f %6.1f %.4e %.4e %.6f %.3f",i, D[ibin], R, z, sigW, sigzadd, fvBH, zd0, zdBH, rhoMS, rhoBH, rhoBH / rhoMS, fBHtmp);
        fBHtmp *= rhoBH / rhoMS;
      }
      fBH[i][ibin] = fBHtmp;  // later multiplied by addwtj
      if (i < 8 && printBHfac){
        printf(" %.6f\n",fBH[i][ibin]);
      }
      double nMS = (i == 8) ? n0MSb*rhos[8] : (i == 9) ? n0MSND*rhos[9] : (i == 10) ? n0MSSH*rhos[10] : n0MSd[i]*rhos[i];
      double rho = (i == 8) ? n0b  *rhos[8] : (i == 9) ? n0ND  *rhos[9] : (i == 10) ? n0SH  *rhos[10] : n0d[i]  *rhos[i];
      if (AI0 > 0 && Isen - Isst > 0 && EVI0 > 0 && VIsen - VIsst > 0){
        double fIVIs = fIVI_detect(extI, Isst, Isen, extVI, VIsst, VIsen, i);
        rhoD_S[i][ibin] = nMS * fIVIs * 1e-06 * D[ibin] * D[ibin];
        nallS += rhoD_S[i][ibin]*dD;
        // printf (" %.5f %.5e",fIVIs,rhoD_S[i][ibin]);
      }else if (AI0 > 0 && Isen - Isst > 0){
        double fIs = fLF_detect(extI, Isst, Isen, i);
        rhoD_S[i][ibin] = nMS * fIs * 1e-06 * D[ibin] * D[ibin];
        nallS += rhoD_S[i][ibin]*dD;
        // printf (" %.5f %.5e",fIs,rhoD_S[i][ibin]);
      }else{
        double tmpDswt = (gammaDs == 0.5) ? sqrt(D[ibin]/8000.0)  // sqrt = 2.0 (volume effect) - 1.5 (limiting mag effect), ideally LF(I) & AI should be used 
                       : pow(((D[ibin]+10)/8000.0), fabs(gammaDs));
        if (gammaDs < 0) tmpDswt = 1 / tmpDswt;
        rhoD_S[i][ibin] = nMS * tmpDswt * 1e-03;
      }
      if (wtD_L != 0) rho *= pow((D[ibin] + 1000)/4500. , wtD_L);
      rhoD_L[i][ibin] = (CheckD == 1) ? rho * 1e-06 * D[ibin] * D[ibin] : rho;
      // Cumulative. The min and max edges are not correctly treated. 
      // But practically OK especially considering a risk of ran < Cumu_min or ran > Cumu_max when strictly considered.
      cumu_rho_S[i][ibin]  = (ibin==0) ? 0 : cumu_rho_S[i][ibin-1] + 0.5*(rhoD_S[i][ibin-1] + rhoD_S[i][ibin])*dD;
      cumu_rho_L[i][ibin]  = (ibin==0) ? 0 : cumu_rho_L[i][ibin-1] + 0.5*(rhoD_L[i][ibin-1] + rhoD_L[i][ibin])*dD;
      cumu_rho_all_S[ibin] += cumu_rho_S[i][ibin];
      cumu_rho_all_L[ibin] += cumu_rho_L[i][ibin];
      rhosum += (printrhoS) ? rhoD_S[i][ibin] : rho;
      if (ibin%npri==0){
        if (printrhoS){
          printf (" %d: %.1e ",i,rhoD_S[i][ibin]);
          printf ("( %.2e )",cumu_rho_S[i][ibin]);
        }else{
          printf (" %d: %.1e ",i,rho);
          printf ("( %.2e )",cumu_rho_L[i][ibin]);
        }
      }
    }
    // printf ("\n");
    if (ibin%npri==0){
      printf (" All: %.1e ",rhosum);
      if (printrhoS){
        printf ("( %.2e )\n",cumu_rho_all_S[ibin]);
      }else{
        printf ("( %.2e )\n",cumu_rho_all_L[ibin]);
      }
    }
  }
  if (printBHfac)
    exit(1);
  int **ibinptiles_S, **ibinptiles_L;
  ibinptiles_S  = (int **)malloc(sizeof(int *) * ncomp);
  ibinptiles_L  = (int **)malloc(sizeof(int *) * ncomp);
  for (int i=0; i<ncomp; i++){
    ibinptiles_S[i] = (int *)calloc(22, sizeof(int *));
    ibinptiles_L[i] = (int *)calloc(22, sizeof(int *));
  }
  for (int i=0;i<ncomp;i++){
    // Store percentiles
    double norm_S = cumu_rho_S[i][nbin];
    double norm_L = cumu_rho_L[i][nbin];
    if (norm_S == 0 && i == 9) continue; // when NSD is not considered
    if (norm_S == 0 && i == 10) continue; // when Stellar halo is not considered
    for (int ibin=0; ibin<=nbin;ibin++){
      double Pnorm_S = cumu_rho_S[i][ibin] / norm_S;
      int intp_S = Pnorm_S*20;
      if (ibinptiles_S[i][intp_S] == 0) ibinptiles_S[i][intp_S] = (intp_S==0) ? 1 : ibin+0.5;
      double Pnorm_L = cumu_rho_L[i][ibin] / norm_L;
      int intp_L = Pnorm_L*20;
      if (ibinptiles_L[i][intp_L] == 0) ibinptiles_L[i][intp_L] = (intp_L==0) ? 1 : ibin+0.5;
    }
  }
  
  // Calculate optical depth when Isen - Isst && AIrc are given (not needed for Monte Carlo simulation)
  int CALCTAU = getOptioni(argc,argv, "CALCTAU", 1, 0);
  double tauall = 0, Nsall = 0;
  if (CALCTAU && Isen - Isst > 0 && AIrc > 0){
    void calc_opticaldepth(double *tauall, double *Nsall, int idata, int Dsmax21, double AI0, double hscale, double Isst, double Isen);
    calc_opticaldepth(&tauall, &Nsall, idata, Dmax, AI0, hscale, Isst, Isen);
  }

  // check D distribution
  if (CheckD == 1){
    double getcumu2xist (int n, double *x, double *F, double *f, double Freq, int ist, int inv);
    for (int i=0; i<500000; i++){
      double ran = ran1();
      double cumu = 0;
      int j_L, j_S;
      for (j_L=0;j_L<ncomp;j_L++){
        cumu += cumu_rho_L[j_L][nbin]/cumu_rho_all_L[nbin];
        if (ran < cumu) break;
      }
      cumu = 0;
      for (j_S=0;j_S<ncomp;j_S++){
        cumu += cumu_rho_S[j_S][nbin]/cumu_rho_all_S[nbin];
        if (ran < cumu) break;
      }
      ran = ran1();
      int inttmp = ran*20;
      int kst = 1;
      for (int itmp = inttmp; itmp > 0; itmp--){
        kst = ibinptiles_L[j_L][itmp];
        if (kst > 0) break;
      }
      double d_L = getcumu2xist(nbin+1, D, cumu_rho_L[j_L],rhoD_L[j_L],ran*cumu_rho_L[j_L][nbin],kst,0);
      kst = 1;
      for (int itmp = inttmp; itmp > 0; itmp--){
        kst = ibinptiles_S[j_S][itmp];
        if (kst > 0) break;
      }
      double d_S = getcumu2xist(nbin+1, D, cumu_rho_S[j_S],rhoD_S[j_S],ran*cumu_rho_S[j_S][nbin],kst,0);
      // printf ("%d %6.0f %d %6.0f\n",j_S,d_S,j_L,d_L);
      // pick velocities
      void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD); //
      double vxyz_L[3] = {};
      double tau_l = (j_L == 9) ? mageND + sageND*gasdev()
                   : (j_L == 8) ? mageB + sageB*gasdev() 
                   : (j_L ==10) ? 14 
                   : medtauds[j_L];
      get_vxyz_ran(vxyz_L, j_L, tau_l, d_L, lDs[idata], bDs[idata]);
      double vx_l = vxyz_L[0];
      double vy_l = vxyz_L[1];
      double vz_l = vxyz_L[2];
      double v_l = sqrt(vx_l*vx_l + vy_l*vy_l + vz_l*vz_l);

      printf ("%d %6.0f %d %6.0f %.6f %7.2f %7.2f %7.2f %7.2f\n",j_S,d_S,j_L,d_L, sqrt(d_L/8000.0), vx_l, vy_l, vz_l, v_l);
    }
    exit(1);
  }

  // release luminosity functions
  population.release_luminosity_functions();

  // Monte Carlo simulation
  printf ("#----- Output of Monte Carlo simulation w/ VERBOSITY= %d and seed= %ld -------- \n",VERBOSITY,seed);
  double getcumu2xist (int n, double *x, double *F, double *f, double Freq, int ist, int inv);
  MonteCarloStats stats;
#define ncntall stats.ncntall
#define ncnts stats.ncnts
#define ncntbWD stats.ncntbWD
#define ncntbCD stats.ncntbCD
#define nBD stats.nBD
#define nMS stats.nMS
#define nWD stats.nWD
#define nNS stats.nNS
#define nBH stats.nBH
#define Nlike stats.Nlike
#define NrejIS stats.NrejIS
#define Ngen stats.Ngen
#define wtlike stats.wtlike
#define wtlike_tE stats.wtlike_tE
#define wtlike_except_piE stats.wtlike_except_piE
#define wtlike_w_piEe stats.wtlike_w_piEe
#define wtlike_except_thE stats.wtlike_except_thE
#define wtlike_w_thEe stats.wtlike_w_thEe
#define SumGamma stats.SumGamma
#define SumtE stats.SumtE
#define logtEmin stats.logtEmin
#define logtEmax stats.logtEmax
#define NbintE stats.NbintE
#define NlogtEs stats.NlogtEs
  if (VERBOSITY == 2) printf ("#        wtj           tE       thetaE          piEN          piEE   D_S         muSl         muSb iS iL fREM");
  if (VERBOSITY == 3) printf ("#        wtj          M_L   D_L   D_S          t_E      theta_E         pi_E         pi_EN         pi_EE       mu_rel        mu_Sl        mu_Sb     I_L     K_L iS iL fREM");
  if (VERBOSITY == 4) printf ("#        wtj          M_L   D_L   D_S          t_E      theta_E         pi_E         pi_EN         pi_EE       mu_rel        mu_Sl        mu_Sb     I_L     K_L iS iL fREM   muhel_N      muhel_E        muhel");
  if (VERBOSITY == 5) printf ("#        wtj          M_L   D_L   D_S          t_E      theta_E         pi_E         pi_EN         pi_EE       mu_rel        mu_Sl        mu_Sb     I_L     K_L iS iL fREM   muhel_N      muhel_E        muhel      R_E");
  if (VERBOSITY >= 2 && BINARY    == 1) printf ("     q21         M2         aL     aLpmin         u0 BL");
  if (VERBOSITY >= 2) printf ("\n");
  for (long j=0; j< NSIMU; j++){
     if (j==NSIMU-1 && NlikeMIN > 0){
       long NSIMUnew = (NSIMU-1) * (double) NlikeMIN/(Nlike+1); // +1 for case of N_like = 0
       if (NSIMUnew > NSIMU) printf ("# Increase NSIMU to %ld\n",NSIMUnew);
       if (NSIMUnew > NSIMU) NSIMU = NSIMUnew;
       if (NSIMU > NDATAMAX) printf ("# Decrease NSIMU to %ld\n",NDATAMAX) ;
       if (NSIMU > NDATAMAX) NSIMU = NDATAMAX;
     }
     double ran, cumu, addGammaIS = 1;
     int inttmp, kst;
     Ngen++;

     // pick D_s
     ran = ran1(); 
     cumu = 0;
     int i_s;
     for (i_s=0;i_s<ncomp;i_s++){
        cumu += cumu_rho_S[i_s][nbin]/cumu_rho_all_S[nbin];
        if (ran < cumu) break;
     }
     if (i_s == ncomp){ // Sometimes happened
       j--;
       continue; 
     }
     double tau_s = (i_s == 9) ? mageND + sageND*gasdev()
                  : (i_s == 8) ? mageB + sageB*gasdev() 
                  : (i_s ==10) ? 14 
                  : medtauds[i_s];
     ran = ran1();
     inttmp = ran*20;
     kst = 1;
     for (int itmp = inttmp; itmp > 0; itmp--){
       kst = ibinptiles_S[i_s][itmp];
       if (kst > 0) break;
     }
     ran = ran* cumu_rho_S[i_s][nbin];
     double D_s = getcumu2xist(nbin+1, D, cumu_rho_S[i_s],rhoD_S[i_s],ran,kst,0);
     int nbinDs = floor(D_s/dD);

     // pick D_l

     /***************************************************************/
     /*** Importance sampling when thetaErange and piErange are given ***/
     // This part is added on July 8, 2022
     // Determine D_l range when thetaErange and piErange are both given
     int nbinDlmin = 0, nbinDlmax = nbinDs;
     if (thetaEmax - thetaEmin > 0 && piEmax - piEmin > 0){
       double Dlmin = 1000 / (piEmax*thetaEmax + 1000/D_s);
       double Dlmax = 1000 / (piEmin*thetaEmin + 1000/D_s);
       nbinDlmin = floor(Dlmin/dD);
       nbinDlmax = floor(Dlmax/dD) + 1;
       if (nbinDlmin < 0) nbinDlmin = 0; 
       if (nbinDlmax > nbinDs) nbinDlmax = nbinDs; 
     }
     /***************************************************************/
     
     ran = ran1(); 
     cumu = 0;
     int i_l;
     for (i_l=0;i_l<ncomp;i_l++){
       cumu += (cumu_rho_L[i_l][nbinDlmax] - cumu_rho_L[i_l][nbinDlmin])/(cumu_rho_all_L[nbinDlmax] - cumu_rho_all_L[nbinDlmin]);
       if (ran < cumu) break;
     }
     if (i_l == ncomp){ // when nbinDs == 0
       j--;
       continue; 
     }
     double tau_l = (i_l == 9) ? mageND + sageND*gasdev()
                  : (i_l == 8) ? mageB + sageB*gasdev()
                  : (i_l == 10) ? 14
                  : medtauds[i_l];
     ran = ran1() * (cumu_rho_L[i_l][nbinDlmax] - cumu_rho_L[i_l][nbinDlmin]) + cumu_rho_L[i_l][nbinDlmin];
     inttmp = ran*20 / cumu_rho_L[i_l][nbin]; // make inttmp < 1
     kst = 1;
     // printf ("inttmp= %2d\n",inttmp);
     for (int itmp = inttmp; itmp > 0; itmp--){
       kst = ibinptiles_L[i_l][itmp];
       if (kst > 0) break;
     }
     double D_l = getcumu2xist(nbin+1, D, cumu_rho_L[i_l],rhoD_L[i_l],ran,kst,0);
     addGammaIS *=  (cumu_rho_all_L[nbinDlmax] - cumu_rho_all_L[nbinDlmin])/ cumu_rho_all_L[nbin];
     // if (i_l == 10)
     //   printf("ran= %.5e inttmp= %2d kst= %d D_l= %f\n",ran,inttmp, kst, D_l);


     // pick velocities
     void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD); //
     double vxyz_S[3] = {}, vxyz_L[3] = {};
     get_vxyz_ran(vxyz_S, i_s, tau_s, D_s, lDs[idata], bDs[idata]);
     get_vxyz_ran(vxyz_L, i_l, tau_l, D_l, lDs[idata], bDs[idata]);
     double vx_s = vxyz_S[0], vx_l = vxyz_L[0];
     double vy_s = vxyz_S[1], vy_l = vxyz_L[1];
     double vz_s = vxyz_S[2], vz_l = vxyz_L[2];

     // Convert into proper motion
     // use relation: (not exactly correct if zsun > 0)
     //   v_l = v_x sinl      + v_y cosl
     //   v_b = v_x cosl sinb - v_y sinl sinb + v_z cosb
     double vxrel_s = vx_s - vxsun;
     double vyrel_s = vy_s - vysun;
     double vzrel_s = vz_s - vzsun;
     double muSl   = (vxrel_s*sinl      + vyrel_s*cosl)*KS2MY/D_s;  // iru
     double muSb   = (vxrel_s*cosl*sinb - vyrel_s*sinl*sinb + vzrel_s*cosb)*KS2MY/D_s; // iru
     double murells[3] = {}, murelbs[3] = {}, murels[3] = {}; // 0: MS/WD, 1: NS, 2: BH
     double thetakick = (REMNANT == 1 && onlyWD == 0) ? asin(1 - 2*ran1()) : 0;
     double phikick   = (REMNANT == 1 && onlyWD == 0) ? ran1()*2*PI - PI : 0;
     int nREM = (REMNANT == 1 && onlyWD == 0) ? 3 : 1;
     for (int iREM = 0; iREM < nREM; iREM++){
       double vxadd = 0, vyadd = 0, vzadd = 0;
       if (iREM > 0){ // iREM == 1 or 2
         double vkick = (iREM == 1) ? vkickNS : vkickBH;  // Table 2 of Lam et al. 2020
         if (MXDkick == 1){ // maxwell distribution
           double sigv1D = vkick * 0.5 * sqrt(0.5*PI); // ~0.6266 vkick
           vxadd = sigv1D * gasdev();
           vyadd = sigv1D * gasdev();
           vzadd = sigv1D * gasdev();
         }else{ // fixed vkick = vavg
           vxadd =  vkick * cos(thetakick) * cos(phikick);
           vyadd =  vkick * cos(thetakick) * sin(phikick);
           vzadd =  vkick * sin(thetakick);
         }
       }
       double vxrel_l = vx_l + vxadd - vxsun;
       double vyrel_l = vy_l + vyadd - vysun;
       double vzrel_l = vz_l + vzadd - vzsun;
       double muLl   = (vxrel_l*sinl      + vyrel_l*cosl)*KS2MY/D_l;
       double muLb   = (vxrel_l*cosl*sinb - vyrel_l*sinl*sinb + vzrel_l*cosb)*KS2MY/D_l;
       double murellhel = muLl - muSl;
       double murelbhel = muLb - muSb;
       murells[iREM] = murellhel - vEarthl*KS2MY*(D_s-D_l)/D_s/D_l; // hel -> geo, iru
       murelbs[iREM] = murelbhel - vEarthb*KS2MY*(D_s-D_l)/D_s/D_l; // hel -> geo, iru
       murels[iREM]  = sqrt(murells[iREM]*murells[iREM] + murelbs[iREM]*murelbs[iREM]);
     }
     double murell = murells[0]; // These are updated when fREM >= 2 later
     double murelb = murelbs[0];
     double murel  = murels[0];
     

     // Lens mass
     /***************************************************************/
     /*** Importance sampling when thetaErange or piErange is given ***/
     // This part is added on July 6, 2022

     // Calculate pirel 
     double pirel = 1000*(1/D_l - 1/D_s);

     // Determine max present day mass
     double MBHmax = 15.6775; // max mass when REMNANT when Mini ~ 42.21, depend on Mini2Mrem!!!
     double MWDmax = 1.375; // max mass for WD
     double Minidie; // max mass for MS
     if (i_l == 8){ // bulge 
       int iage_l = tau_l * 2 + 0.5;
       iage_l *= 50;
       int itmp = (iage_l - agesB[0])/(agesB[1] - agesB[0]);
       Minidie = MinidieB[itmp];
     }else if(i_l == 10){ // Halo
       Minidie = MinidieD[nageD-1]; // nageD-1: halo
     }else if(i_l == 9){ // NSD
       int iage_l = 100*(tau_l + 0.5);
       int itmp = (iage_l <= agesND[0]) ? 0
                : (iage_l >= agesND[nageND-1]) ? nageND-1 
                : (iage_l - agesND[0])/(agesND[1] - agesND[0]);
       Minidie = MinidieND[itmp];
     }else if(i_l == 7){ // thick disk
       Minidie = MinidieD[nageD-2]; // nageD-1: halo
     }else{ // thin disk
       int iage_l = tau_l * 100 + 0.5;
       iage_l = (iage_l % 5 > 2.5) ? iage_l + (5 - iage_l % 5) : iage_l - iage_l % 5;
       if (iage_l < 5) iage_l = 5;
       int itmp = (iage_l - agesD[0])/(agesD[1] - agesD[0]);
       Minidie = MinidieD[itmp];
     }
     double MMSmax = (BINARY == 1) ? 2 * Minidie : Minidie; // equal-mass binary is max mass when BINARY
     double MPDmax = (REMNANT == 1) ? MBHmax 
                   : (onlyWD == 1 && MWDmax > MMSmax) ? MWDmax
                   : MMSmax;
     double MPDmin = Ml;
     
     // Calculate Mrange from tErange, thetaErange and piErange
     double MtEmin = MPDmax, MtEmax = MPDmin, MthEmin = MPDmin, MthEmax = MPDmax, MpiEmin = MPDmin, MpiEmax = MPDmax;
     if (tEmax - tEmin > 0){
       double MMSWDmax = (MWDmax > MMSmax && (onlyWD == 1 || REMNANT == 1)) ? MWDmax : MMSmax;
       double tmpfac = murel * murel / KAPPA/pirel/365.25/365.25;
       double MtEMSmin = tmpfac*tEmin*tEmin;
       if (MtEMSmin < MPDmin) MtEMSmin = MPDmin;
       double MtEMSmax = tmpfac*tEmax*tEmax;
       if (MtEMSmax > MMSWDmax) MtEMSmax = MMSWDmax;
       if (MtEMSmax > MPDmin && MtEMSmin < MMSWDmax){
         MtEmin = MtEMSmin;
         MtEmax = MtEMSmax;
       }
       if (REMNANT == 1 && onlyWD == 0){
         // With NS kick velocity
         tmpfac = murels[1] * murels[1] / KAPPA/pirel/365.25/365.25; // NS
         double MtENSmin = tmpfac*tEmin*tEmin;
         if (MtENSmin < MNSMIN) MtENSmin = MNSMIN;
         double MtENSmax = tmpfac*tEmax*tEmax;
         if (MtENSmax > MNSMAX) MtENSmax = MNSMAX;
         if (MtENSmax > MNSMIN && MtENSmin < MNSMAX){
           MtEmin = (MtENSmin < MtEmin) ? MtENSmin : MtEmin;
           MtEmax = (MtENSmax > MtEmax) ? MtENSmax : MtEmax;
         }
         // With BH kick velocity
         double MBHmin = 4.9911; // when Mini ~ 15 Msun
         tmpfac = murels[2] * murels[2] / KAPPA/pirel/365.25/365.25; // BH
         double MtEBHmin = tmpfac*tEmin*tEmin;
         if (MtEBHmin < MBHmin) MtEBHmin = MBHmin;
         double MtEBHmax = tmpfac*tEmax*tEmax;
         if (MtEBHmax > MBHmax) MtEBHmax = MBHmax;
         if (MtEBHmax > MBHmin && MtEBHmin < MBHmax){
           MtEmin = (MtEBHmin < MtEmin) ? MtEBHmin : MtEmin;
           MtEmax = (MtEBHmax > MtEmax) ? MtEBHmax : MtEmax;
         }
         // printf("# MtEMSmin-max= %.6f - %.6f %.6f\n",MtEMSmin,MtEMSmax,murels[0]);
         // printf("# MtENSmin-max= %.6f - %.6f %.6f\n",MtENSmin,MtENSmax, murels[1]);
         // printf("# MtEBHmin-max= %.6f - %.6f %.6f\n",MtEBHmin,MtEBHmax, murels[2]);
         // printf("# MtEmin-max  = %.6f - %.6f\n",MtEmin,MtEmax);
       }
       if (MtEmin > MtEmax){
         j--;
         NrejIS++;
         continue;
       }
     }else{
       MtEmin = MPDmin; // The initial values were opposite to define MtEmin and MtEmax when tErange is given
       MtEmax = MPDmax;
     }
     if (thetaEmax - thetaEmin > 0){
       double Mtmp;
       Mtmp = thetaEmin*thetaEmin/KAPPA/pirel;
       MthEmin = (Mtmp < MPDmin) ? MPDmin 
               : (Mtmp > MPDmax) ? MPDmax
               : Mtmp;
       Mtmp = thetaEmax*thetaEmax/KAPPA/pirel;
       MthEmax = (Mtmp < MPDmin) ? MPDmin 
               : (Mtmp > MPDmax) ? MPDmax
               : Mtmp;
     }
     if (piEmax - piEmin > 0){
       double Mtmp;
       Mtmp = pirel/KAPPA/piEmax/piEmax; // lower value
       MpiEmin = (Mtmp < MPDmin) ? MPDmin 
               : (Mtmp > MPDmax) ? MPDmax
               : Mtmp;
       Mtmp = pirel/KAPPA/piEmin/piEmin; // larger value
       MpiEmax = (Mtmp < MPDmin) ? MPDmin 
               : (Mtmp > MPDmax) ? MPDmax
               : Mtmp;
     }

     // Reject when no overlap between allowed Mranges by thetaE and piE
     if (MtEmin >= MthEmax || MtEmin >= MpiEmax || MthEmin >= MpiEmax || 
         MtEmax <= MthEmin || MtEmax <= MpiEmin || MthEmax <= MpiEmin || 
         MtEmax - MtEmin <= 0 || MthEmax - MthEmin <= 0 || MpiEmax - MpiEmin <= 0)
     {
       j--;
       NrejIS++;
       continue;
     }

     // Determine Mrange
     int nbinMmin = 0, nbinMmax = 0;
     if (tEmax - tEmin > 0 || thetaEmax - thetaEmin > 0 || piEmax - piEmin > 0){
       double Minimin = (MtEmin > MthEmin) ? MtEmin : MthEmin;
       if (MpiEmin > Minimin) Minimin = MpiEmin;
       double Minimax = (MtEmax < MthEmax) ? MtEmax : MthEmax;
       if (MpiEmax < Minimax) Minimax = MpiEmax;
       if (REMNANT == 1 && onlyWD == 0){ // When can be NS or BH
         double MWDmin = 0.109 * Minidie + 0.394;
         // max: 
         // 1. If NS is possible, consider up to Mu (120 Msun) cuz Mrem2Mini is difficult for NS or BH.
         // 2. If NS is impossible but WD is possible, consider up to the corresponding initial mass.
         // 3. If NS or WD is impossible (i.e., star), initial mass = present-day mass
         Minimax = (Minimax > MNSMIN) ? Mu  // Consider up to Mu when exceed min mass of NS (~)
                 : (Minimax > MWDmin) ? (Minimax - 0.394)/0.109  // 
                 : Minimax;
         // min: 
         // 1. If need to be NS/BH, minimum mass is MiniWDmax
         // 2. If need to be WD, minimum mass is the corresponding initial mass.
         // 3. If not need to be WD/NS/BH, initial mass = present-day mass 
         Minimin = (Minimin > MWDmax && Minimin > MMSmax) ? MiniWDmax // need to be NS or BH
                 : (Minimin > MMSmax) ? (Minimin - 0.394)/0.109  // Minidie < Minimin < MWDmax, need to be WD
                 : (BINARY == 1) ? 0.5 * Minimin // equal mass binary
                 : Minimin;
       }else if(onlyWD == 1){  // When can be WD
         double MWDmin = 0.109 * Minidie + 0.394;
         Minimax = (Minimax > MWDmin) ? (Minimax - 0.394)/0.109  // can be WD
                 : Minimax;
         Minimin = (Minimin > MMSmax) ? (Minimin - 0.394)/0.109  // need to be WD
                 : (BINARY == 1) ? 0.5 * Minimin // equal mass binary
                 : Minimin;
       }else{
         Minimin = (BINARY == 1) ? 0.5 * Minimin // equal mass binary
                 : Minimin;
       }
       nbinMmin = floor((log10(Minimin) - logMst)/dlogM);
       nbinMmax = floor((log10(Minimax) - logMst)/dlogM) + 1;
       if (nbinMmin < 0) nbinMmin = 0;
       if (nbinMmax > nm) nbinMmax = nm;
       // printf("# MtEmin-max = %.6f - %.6f\n",MtEmin,MtEmax);
       // printf("# MthEmin-max= %.6f - %.6f\n",MthEmin,MthEmax);
       // printf("# MpiEmin-max= %.6f - %.6f\n",MpiEmin,MpiEmax);
       // printf("# Minimin-max= %.6f - %.6f\n",Minimin,Minimax);
     }

     /***************************************************************/
     
     double logM, M_l;
     if (nbinMmax - nbinMmin > 0){  // importance sampling
       ran = ran1() * (PlogM_cum_norm_B[nbinMmax] - PlogM_cum_norm_B[nbinMmin]) + PlogM_cum_norm_B[nbinMmin]; // limit Mmin - Mmax
       addGammaIS *= PlogM_cum_norm_B[nbinMmax] - PlogM_cum_norm_B[nbinMmin]; // reduce weight by the confined region
     }else{
       ran = ran1();
     }
     inttmp = ran*20;
     kst = 1; // to avoid bug when inttmp = 0
     for (int itmp = inttmp; itmp > 0; itmp--){
       kst = imptiles_B[itmp] - 1; 
       if (kst > 0) break;
     }
     logM = getcumu2xist(nm, logMass_B, PlogM_cum_norm_B, PlogM_B, ran, kst, 0);
     M_l = pow(10, logM);
     double Morg = M_l; // for wtM_L

     // Reject or Evolve into WD
     int fREM = 0;
     void Mini2Mrem (double *pout, double M, int mean); 
     // printf "# iage_l M_l > Minidie{iage_l}" if M_l > Minidie{iage_l}; 
     if (M_l > Minidie){
       if (REMNANT==1 || onlyWD == 1){
         double pout[2] = {};
         Mini2Mrem(pout, M_l, 0);  // 0 : random
         M_l  = pout[0]; // Mass after evolution
         fREM = pout[1]; // fREM should be double if mean == 1, but int here cuz mean == 0
         if (fREM >= 2){ // Add kick velocity for NS or BH
            if (onlyWD == 1){ // reject
              j--;
              continue;
            }
            murell = murells[fREM - 1];
            murelb = murelbs[fREM - 1];
            murel  = murels[fREM - 1];
         }
       }else{ // reject
         j--;
         continue;
       }
     }

     // Binary system assuming the picked M is a primary
     // -- Binary distribution developed by Koshimoto+2020, AJ, 159, 268 is used
     int swl = 0; // 0: single, 1: close binary, 2: wide binary
     double M_l2 = 99999, q21 = 99999, u0S = 99999, al = 99999, alpmin = 99999, apdetS = 99999, apdetL = 99999;
     if (BINARY && fREM == 0){  // Remnant in a binary should be ideally considered, but currently not yet
       double mult = 0.196 + 0.255*M_l; // Table 2 of Koshimoto+20, AJ, 159, 268
       if (mult > MAXMULT) mult = MAXMULT;
       ran = ran1();
       double coeff, q2;
       if (ran < 0.5 * mult){ // close binary
         swl = 1; // 0: single, 1: close binary, 2: wide binary
         double gamma =  1.16 - 2.79*log10(M_l); // Table 2 of Koshimoto+20, AJ, 159, 268
         if (gamma > MAXGAMMA) gamma = MAXGAMMA;
         if (gamma < MINGAMMA) gamma = MINGAMMA;
         coeff = -1;
         double tmp = pow(0.1, gamma+1); // because we ignore q < 0.1
         q2 = pow( (1-tmp)*ran1() + tmp, 1/(gamma+1) ); // inverse transform sampling
       }else if(ran < mult){
         swl = 2; // 0: single, 1: close binary, 2: wide binary
         double gamma = (M_l>=0.344) ? 0 : -3.09 - 6.67*log10(M_l); // Table 2 of Koshimoto+20, AJ, 159, 268 
         if (gamma > MAXGAMMA) gamma = MAXGAMMA;
         if (gamma < MINGAMMA) gamma = MINGAMMA;
         coeff = 1;
         double tmp = pow(0.1, gamma+1); // because we ignore q < 0.1
         q2 = pow( (1-tmp)*ran1() + tmp, 1/(gamma+1) ); // inverse transform sampling
       }
       if (swl > 0){
         M_l2 = M_l * q2;
         // Determine which is the main lens (i.e., arbitrary star in Koshimoto+2020)
         double Ptmp = sqrt(M_l)/(sqrt(M_l) + sqrt(M_l2)); // ratio of lensing zone size
         ran = ran1();
         double Mtmp1 = (ran < Ptmp) ? M_l  : M_l2; // lens candidate
         double Mtmp2 = (ran < Ptmp) ? M_l2 : M_l;  // lens companion candidate, could be lens together
         double qtmp = Mtmp2/Mtmp1; // Can be > 1
         double thetaE1  = sqrt(Mtmp1* pirel*KAPPA);
         double thetaE12 = sqrt((Mtmp1+Mtmp2)* pirel*KAPPA);
         //   pick up aproj
         double pout[2] = {};
         void getaproj(double *pout, double M1, double M2, int coeff);
         getaproj(pout, M_l, M_l2, coeff);
         al     = (pout[0] < 99) ? pow(10.0, pout[0]) : -1;
         alpmin = pout[1];
         //   set upper and lower limit for sep based on u0 and 
         //   central caustic size (Chung+2005, ApJ, 630, 535)
         u0S = (u0obs > 0) ? u0obs : ran1();  // u0 following uniform distribution if not given
         double tmp = qtmp/u0S; 
         apdetL = (sqrt(tmp + 1) + sqrt(tmp))*thetaE1*D_l*0.001; //w < u0S using w=4q/(s-1/s)^2
         if (qtmp > 1) tmp = 1/qtmp/u0S; // because q cannot be over 1 when lens is close binary
         apdetS = (sqrt(tmp + 1) - sqrt(tmp))*thetaE12*D_l*0.001; //w < u0S using w=4q/(s-1/s)^2
         apdetS *= (qtmp > 1) ? sqrt(1+1/qtmp) : sqrt(1+qtmp); // 
         if (alpmin > apdetL){  // Lens = arbitrary star but with undetectable wide companion
           swl += 10; // Undetectable 
           M_l  = Mtmp1;
           M_l2 = Mtmp2;
           q21 = qtmp;
         }else if (alpmin < apdetS){ // Lens = very close binary which looks like a single lens
           swl += 20; // Undetectable 
           M_l  = Mtmp1 + Mtmp2;
           M_l2 = 99;
           q21 = qtmp;
         }else{ // Detectable binary-lens event
           j--;
           continue;
         }
       }
     }

     // Microlens parameters
     double thetaE = sqrt(M_l* pirel * KAPPA); // 99 just for DSGAMMA
     if (thetaE == 0){ // oukyuu sochi
       printf ("#ERROR: thetaE = 0!! M_l=  %.5e\n",M_l);
       j--;
       continue;
     }
     double tE = thetaE/murel*365.25;
     double piE  = pirel/thetaE;
     double piEN = piE * ( murelb*cosPA + murell*sinPA)/murel;
     double piEE = piE * (-murelb*sinPA + murell*cosPA)/murel;
     double Gamma = 1.6e-08 * D_l*D_l*thetaE*murel; // 1.6e-08 makes Gamma to be < ~1
     if (NoGAMMAIS == 1){
        Gamma  = 8e-09 * D_l*D_l*thetaE*murel; // 8e-09 makes Gamma to be < ~1 when addGammaIS multiplied
        Gamma *= addGammaIS;
        addGammaIS = 1;
     }
     // printf("Ml= %.6f, D_l= %.2f Gamma= %.5e piE= %.5f thetaE= %.5f tE= %.3f\n",M_l, D_l,Gamma,piE,thetaE,tE);

     // correction of Gamma for BH
     int ien = (BHhb == 1) ? 9 : 8;
     if (BHhd == 1 && i_l < ien && fREM == 3){
        double fcorBH = interp_x(nbin+1, fBH[i_l], 0, dD, D_l);
        Gamma *= fcorBH;
        // printf("%d %.0f %.6f %.3f %.5e\n",i_l,D_l,fcorBH, murel, Gamma);
     }

     // for mean tE
     SumGamma += Gamma*addGammaIS;
     SumtE    += Gamma*addGammaIS*tE;
     // for median tE
     int ilogtE = (log10(tE) - logtEmin)/(logtEmax - logtEmin)*NbintE;
     if (ilogtE < 0) ilogtE = 0;
     if (ilogtE > NbintE - 1) ilogtE = NbintE - 1;
     NlogtEs[ilogtE] += Gamma*addGammaIS;
     // if (i_l == 10){
     //   printf("j= %ld i_l: %d D_l= %f thetaE= %f murel= %f Gamma= %12.6e\n",j,i_l,D_l, thetaE, murel, Gamma);
     // }
     if (Gamma < ran1() && SMALLGAMMA == 0){
       j--;
       continue;
     }
     // if (i_l == 10)
     //   printf("2 j= %ld i_l: %d D_l= %f thetaE= %f murel= %f Gamma= %12.6e\n",j,i_l,D_l, thetaE, murel, Gamma);
     double like_obs(double mod, double obs, double err, double fe, int det, int uniform_mode);
     double addGamma = 1;
     double wtj = (Gamma > 1 || SMALLGAMMA == 1) ? Gamma : 1;

     // Correction for wtM_L and wtD_L
     if (wtM_L != 0)  addGamma /= pow(Morg/0.1 , wtM_L);
     if (wtD_L != 0)  addGamma /= pow((D_l + 1000)/4500. , wtD_L);

     // Constraint from tE
     double Gamma_tE = like_obs(tE, tEobs, tEe, fetE, tEdet, UNIFORM);
     int like_tE = (Gamma_tE > 0) ? 1 : 0;
     if (Gamma_tE > 0) addGamma *= Gamma_tE;
 
     // Constraint from thetaE
     double Gamma_thetaE = like_obs(thetaE, thetaEobs, thetaEe, fethetaE, thetaEdet, UNIFORM);
     int like_thetaE = (Gamma_thetaE > 0) ? 1 : 0;
     if (Gamma_thetaE > 0) addGamma *= Gamma_thetaE;

     // Constraint from piE
     int like_piE = 1;
     if(piEe >0){
       double Gamma_piE = like_obs(piE, piEobs, piEe, fepiE, piEdet, UNIFORM);
       like_piE = (Gamma_piE > 0) ? 1 : 0;
       if (Gamma_piE > 0) addGamma *= Gamma_piE;
     }else if(piENe > 0 && piEEe > 0){
       double Gamma_piEN = like_obs(piEN, piENobs, piENe, fepiEN, 0, UNIFORM);
       int like_piEN = (Gamma_piEN > 0) ? 1 : 0;
       if (Gamma_piEN > 0) addGamma *= Gamma_piEN;
       double Gamma_piEE = like_obs(piEE, piEEobs, piEEe, fepiEE, 0, UNIFORM);
       int like_piEE = (Gamma_piEE > 0) ? 1 : 0;
       if (Gamma_piEE > 0) addGamma *= Gamma_piEE;
       like_piE = like_piEN * like_piEE;
     }

     // Constraint from mus (measured in heliocentric)
     int like_mus = 1;
     if (musRCG == 1){ // Recalculate muS relative to RCG
        double musfac  = 1.0 - D_s/Dmean; // need to subtract proper motion of RCG if the origin is RCG. assume v_RCG=0 and D_RCG=Dmean
        vxrel_s = vx_s - musfac*vxsun;
        vyrel_s = vy_s - musfac*vysun;
        vzrel_s = vz_s - musfac*vzsun;
        muSl   = (vxrel_s*sinl      + vyrel_s*cosl)*KS2MY/D_s;
        muSb   = (vxrel_s*cosl*sinb - vyrel_s*sinl*sinb + vzrel_s*cosb)*KS2MY/D_s;
     }
     if(musle > 0 && musbe > 0){
       double Gamma_musl = like_obs(muSl, muslobs, musle, femusl, 0, UNIFORM);
       int like_musl = (Gamma_musl > 0) ? 1 : 0;
       if (Gamma_musl > 0) addGamma *= Gamma_musl;
       double Gamma_musb = like_obs(muSb, musbobs, musbe, femusb, 0, UNIFORM);
       int like_musb = (Gamma_musb > 0) ? 1 : 0;
       if (Gamma_musb > 0) addGamma *= Gamma_musb;
       like_mus = like_musl * like_musb;
     }else if(musNe > 0 && musEe > 0){
       double muSN =  muSb*cosPA + muSl*sinPA;
       double muSE = -muSb*sinPA + muSl*cosPA;
       double Gamma_musN = like_obs(muSN, musNobs, musNe, femusN, 0, UNIFORM);
       int like_musN = (Gamma_musN > 0) ? 1 : 0;
       if (Gamma_musN > 0) addGamma *= Gamma_musN;
       double Gamma_musE = like_obs(muSE, musEobs, musEe, femusE, 0, UNIFORM);
       int like_musE = (Gamma_musE > 0) ? 1 : 0;
       if (Gamma_musE > 0) addGamma *= Gamma_musE;
       like_mus = like_musN * like_musE;
     }

     // Constraint from muhel (measured in heliocentric)
     int like_muhel = 1;
     double muhelN, muhelE;
     if ((muhelNe > 0 && muhelEe > 0) || VERBOSITY >= 4){
       double murellhel =  murell + vEarthl*KS2MY*(D_s-D_l)/D_s/D_l; // geo -> hel
       double murelbhel =  murelb + vEarthb*KS2MY*(D_s-D_l)/D_s/D_l; // geo -> hel
       muhelN =  murelbhel*cosPA + murellhel*sinPA;
       muhelE = -murelbhel*sinPA + murellhel*cosPA;
     }
     if(muhelNe > 0 && muhelEe > 0){
       double Gamma_muhelN = like_obs(muhelN, muhelNobs, muhelNe, femuhelN, 0, UNIFORM);
       int like_muhelN = (Gamma_muhelN > 0) ? 1 : 0;
       if (Gamma_muhelN > 0) addGamma *= Gamma_muhelN;
       double Gamma_muhelE = like_obs(muhelE, muhelEobs, muhelEe, femuhelE, 0, UNIFORM);
       int like_muhelE = (Gamma_muhelE > 0) ? 1 : 0;
       if (Gamma_muhelE > 0) addGamma *= Gamma_muhelE;
       like_muhel = like_muhelN * like_muhelE;
     }

     // Constraint from lens brigtness (currently only I and K)
     double getx2y(int n, double *x, double *y, double xin);
     int like_IL = 1;
     double IL = 99;
     if (AI0 > 0 && fREM == 0){
       double M1 = (swl < 20) ? M_l : M_l/(1 + q21);
       IL = (M1 <= 0.010) ? 27.060   //  faintest I mag in input_files/MLemp.dat
          : (M1 >= 1.440) ?  2.355   // brightest I mag in input_files/MLemp.dat
          : getx2y(nMLemp, M_emps, Mag_emps[1], M1);  // 1: I
       if (swl >= 20){
         double M2 = q21*M_l/(1 + q21);
         double IL2 = (M2 <= 0.010) ? 27.060   //  faintest I mag in input_files/MLemp.dat
                    : (M2 >= 1.440) ?  2.355   // brightest I mag in input_files/MLemp.dat
                    : getx2y(nMLemp, M_emps, Mag_emps[1], M2);  // 1: I
         double F12 = pow(10, -0.4*IL) + pow(10, -0.4*IL2);
         IL = -2.5 * log10(F12);
       }
       IL += extinction.distance_modulus_term(D_l) + extinction.at_distance(D_l).i_band; // add DM and AI
       if(ILe > 0){  // This is true when default to cut too bright lens
         double Gamma_IL = like_obs(IL, ILobs, ILe, feIL, ILdet, UNIFORM);
         like_IL = (Gamma_IL > 0) ? 1 : 0;
         if (Gamma_IL > 0) addGamma *= Gamma_IL;
       }
     }

     int like_KL = 1;
     double KL = 99;
     if (AK0 > 0 && fREM == 0){
       double M1 = (swl < 20) ? M_l : M_l/(1 + q21);
       KL = (M1 <= 0.010) ? 30.466   //  faintest K mag in input_files/MLemp.dat
          : (M1 >= 1.440) ?  1.535   // brightest K mag in input_files/MLemp.dat
          : getx2y(nMLemp, M_emps, Mag_emps[4], M1);  // 4: K
       if (swl >= 20){
         double M2 = q21*M_l/(1 + q21);
         double KL2 = (M2 <= 0.010) ? 30.466   //  faintest K mag in input_files/MLemp.dat
                    : (M2 >= 1.440) ?  1.535   // brightest K mag in input_files/MLemp.dat
                    : getx2y(nMLemp, M_emps, Mag_emps[4], M2);  // 4: K
         double F12 = pow(10, -0.4*KL) + pow(10, -0.4*KL2);
         KL = -2.5 * log10(F12);
       }
       KL += extinction.distance_modulus_term(D_l) + extinction.at_distance(D_l).k_band; // add DM and AK
       if(KLe > 0){
         double Gamma_KL = like_obs(KL, KLobs, KLe, feKL, KLdet, UNIFORM);
         like_KL = (Gamma_KL > 0) ? 1 : 0;
         if (Gamma_KL > 0) addGamma *= Gamma_KL;
       }
     }

     // Combine all constraints
     Gamma *= addGamma * addGammaIS;
     wtj   *= addGamma * addGammaIS;
     int like = like_tE*like_thetaE*like_piE*like_mus*like_muhel*like_IL*like_KL;
     int like_expiE = like_tE*like_thetaE*like_mus*like_muhel*like_IL*like_KL;
     int like_exthE = like_tE*like_piE*like_mus*like_muhel*like_IL*like_KL;
     // if (VERBOSITY == 1) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %2d %2d %d %d\n", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, i_s,i_l,fREM,like);
     if (CALCPRIORpiE && piENe > 0 && piEEe > 0 && like_expiE == 1){
       double chi2piEN = (piEN - piENobs)/piENe;
       double chi2piEE = (piEE - piEEobs)/piEEe;
       chi2piEN *= chi2piEN;
       chi2piEE *= chi2piEE;
       wtlike_except_piE += wtj; // bumbo
       wtlike_w_piEe     += wtj * 0.5/PI/piENe/piEEe * exp(-0.5*(chi2piEN + chi2piEE)); // bunshi, sqrt(2*PI*sig^2)*sqrt(2*PI*sig^2)
     } 
     if (CALCPRIORthE && thetaEe > 0 && like_exthE == 1){
       double chi2thE = (thetaE - thetaEobs)/thetaEe;
       chi2thE *= chi2thE;
       wtlike_except_thE += wtj; // bumbo
       wtlike_w_thEe     += wtj * sqrt(0.5/PI)/thetaEe * exp(-0.5*chi2thE); // bunshi, sqrt(2*PI*sig^2)*sqrt(2*PI*sig^2)
     }
     if (like_tE == 1) wtlike_tE += wtj;
     if (like    == 1) Nlike ++, wtlike    += wtj;
     else continue;
     if (VERBOSITY == 1) printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %1d %1d %5.2f %5.2f %d %.5e", wtj,tE,thetaE,piE,M_l,D_s, D_l,vx_s,vy_s,vz_s,vx_l,vy_l,vz_l,murel,i_s,i_l,tau_s,tau_l,fREM,Gamma);
     if (VERBOSITY == 2) printf("%12.6e %12.6e %12.6e %13.6e %13.6e %5.0f %12.5e %12.5e %2d %2d %d", wtj, tE,thetaE,piEN,piEE,D_s,muSl,muSb,i_s,i_l,fREM);
     if (VERBOSITY == 3) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, KL, i_s,i_l,fREM);
     if (VERBOSITY == 4) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d %12.5e %12.5e %12.6e", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, KL, i_s,i_l,fREM, muhelN, muhelE, sqrt(muhelN*muhelN+muhelE*muhelE));
     if (VERBOSITY == 5) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d %12.5e %12.5e %12.6e %.6e", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, KL, i_s,i_l,fREM, muhelN, muhelE, sqrt(muhelN*muhelN+muhelE*muhelE),thetaE*D_l*1e-03);
     if (BINARY == 1 && VERBOSITY > 0)
        printf(" %.4e %.4e %.4e %.4e %.4e %2d",q21, M_l2, al, alpmin, u0S, swl);
     if (VERBOSITY > 0)
        printf ("\n");
     // Count Binary
     ncntall += wtj;
     if (swl == 0) ncnts  += wtj;
     if (swl > 10 && swl < 20) ncntbWD += wtj;
     if (swl > 20) ncntbCD += wtj;
     if (apdetS > apdetL) printf ("# ERROR: apdetS (= %.6f) > apdetL (= %.6f) !!!!!\n",apdetS,apdetL);
     // Count Remnant ()
     if (fREM == 0 && M_l < 0.08) nBD += wtj; // missing BD binaries where M_l (total mass) > 0.08
     if (fREM == 0 && M_l > 0.08) nMS += wtj;
     if (fREM == 1) nWD += wtj;
     if (fREM == 2) nNS += wtj;
     if (fREM == 3) nBH += wtj;
  }
  printf ("# Source number density= %.5e ( %.5e ) arcmin^-2\n",Nsall,nallS*STR2MIN2*1e+6);
  // Median tE (C = 1/2)
  double CumuNlogtE = 0, medtE;
  for (int ilogtE = 0; ilogtE < NbintE; ilogtE++){
    double dP = (ilogtE == 0) ? 0.5*NlogtEs[ilogtE]/SumGamma
              : 0.5*(NlogtEs[ilogtE-1] + NlogtEs[ilogtE])/SumGamma;
    CumuNlogtE += dP;
    if (CumuNlogtE > 0.5){
      double p2 = CumuNlogtE, p1 = CumuNlogtE - dP;
      double dlogtE = (logtEmax - logtEmin)/NbintE;
      double medlogtE = logtEmin + (ilogtE - 0.5)*dlogtE + (0.5 - p1)/(p2 - p1)*dlogtE;
      medtE = pow(10.0, medlogtE);
      // printf("p2= %.5e, p1= %.5e, medlogtE= %.6f, medtE= %.6f\n",p2,p1,medlogtE,medtE);
      break;
    }
  }
  double avetE   = SumtE / SumGamma;
  double everate  = 2 * tauall / PI / avetE * 365.25; // event rate per source per yr
  double everatedeg2 = everate * Nsall * 3600; // event rate per deg^2 per yr (Nsall is in arcmin^2)
  printf ("# avetE= %6.3f days, medtE= %6.3f days, tau= %.6e , event_rate= %.6e /star/yr or %.6e /deg^2/yr\n",avetE,medtE,tauall,everate,everatedeg2);
  
  if (BINARY == 1) printf ("# (n_single n_binwide n_binclose)/n_all= ( %6.0f %6.0f %6.0f ) / %6.0f = ( %.6f %.6f %.6f )\n", ncnts, ncntbWD, ncntbCD, ncntall,ncnts/ncntall,ncntbWD/ncntall,ncntbCD/ncntall);
  printf ("# (n_BD n_MS n_WD n_NS n_BH)/n_all= ( %6.0f %6.0f %6.0f %6.0f %6.0f ) / %6.0f = ( %.6f %.6f %.6f %.6f %.6f )\n", nBD, nMS, nWD, nNS, nBH,ncntall, nBD/ncntall, nMS/ncntall, nWD/ncntall, nNS/ncntall, nBH/ncntall);
  double f_like_tE= wtlike/wtlike_tE;
  if (NrejIS > 0)
    printf ("# N_IS/Ngen= %ld / %ld = %f of events are rejected by importance sampling part\n",NrejIS,Ngen, (double) NrejIS/Ngen);
  printf ("# Nlike/N/Ngen= %d / %ld / %ld     wtlike/wtlike_tE= %.0f / %.0f = %f\n",Nlike,NSIMU,Ngen,wtlike,wtlike_tE,f_like_tE);
  if (CALCPRIORthE)
    printf ("# P_pri(thE)= wtlike_w_thEe/wtlike_except_thE= %.2f / %.2f = %.5e\n",wtlike_w_thEe,wtlike_except_thE,wtlike_w_thEe/wtlike_except_thE);
  if (CALCPRIORpiE)
    printf ("# P_pri(piEN,piEE)= wtlike_w_piEe/wtlike_except_piE= %.2f / %.2f = %.5e\n",wtlike_w_piEe,wtlike_except_piE,wtlike_w_piEe/wtlike_except_piE);
  // gsl_rng_free(r);
  free (D);  
  free (cumu_rho_all_S);
  free (cumu_rho_all_L);
  for (int i=0; i<ncomp; i++){
    free (rhoD_S[i]    );
    free (rhoD_L[i]    );
    free (cumu_rho_S[i]);
    free (cumu_rho_L[i]);
    free (ibinptiles_S[i]);
    free (ibinptiles_L[i]);
  }
  free (rhoD_S    );
  free (rhoD_L    );
  free (cumu_rho_S);
  free (cumu_rho_L);
  free (ibinptiles_S);
  free (ibinptiles_L);
  if (ND == 3){
    for (int i=0; i<nzND; i++){
      for (int j=0; j<nRND; j++){
        free(logsigvNDs[i][j]);
      }
      free(logrhoNDs[i]);
      free(vphiNDs[i]);
      free(corRzNDs[i]);
      free(logsigvNDs[i]);
    }
    free(logrhoNDs);
    free(vphiNDs);
    free(corRzNDs);
    free(logsigvNDs);
  }
  population.release_all();
  free(lDs);
  free(bDs);
  kinematic_tables.release_all();
#undef ncntall
#undef ncnts
#undef ncntbWD
#undef ncntbCD
#undef nBD
#undef nMS
#undef nWD
#undef nNS
#undef nBH
#undef Nlike
#undef NrejIS
#undef Ngen
#undef wtlike
#undef wtlike_tE
#undef wtlike_except_piE
#undef wtlike_w_piEe
#undef wtlike_except_thE
#undef wtlike_w_thEe
#undef SumGamma
#undef SumtE
#undef logtEmin
#undef logtEmax
#undef NbintE
#undef NlogtEs
#undef NSIMU
#undef NlikeMIN
#undef vEarthl
#undef vEarthb
#undef gammaDs
#undef wtD_L
#undef wtM_L
#undef NoGAMMAIS
#undef SMALLGAMMA
#undef VERBOSITY
#undef UNIFORM
#undef BINARY
#undef REMNANT
#undef onlyWD
#undef CALCPRIORpiE
#undef CALCPRIORthE
  return 0;
} // end main

} // namespace genulens
