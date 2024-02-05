/* Generate microlensing events following the Galactic model developed by Koshimoto, Baba & Bennett (2021).
 * N. Koshimoto wrote the original .c version and C. Ranc converted it into .cpp to replace functions (ran1 and gasdev) from the Numerical Recipes in C (NR) with public alternatives (from GSL).
 * We found that the version with ran1 and gasdev from the NR was faster (~1.3 times) than the current version.
 * Please replace them by yourself if you want because we are not allowed to include a NR function in a public package. 
 * Note that a negative seed value has to be used in the NR function in contrast to a positive seed value required for this version. 
 *
 * This is version 1.2 of genulens.
 *
 * Version 1.2 update (Jun 14 - July 11 2022 by N.Koshimoto)
 *   1. NSD (nuclear stellar disk) component is added based on Sormani+22.
 *   2. GC is now located at the position of SgrA* and CenSgrA option is added to locate GC at (l,b) = (0,0)
 *   3. tE mean and median are calculated during the simulation (a large Nsimu value is needed to have a precision).
 *   4. CALCTAU option is added to calculate the optical depth and event rate. (event rate uses the tE mean value based on the simulation.)
 *   5. calc_PA function is added and PA is now automatically calculated from (l, b)
 *   6. Importance sampling used when tE, thetaE, and/or piE is given, which makes most calculation ~10 times faster.
 *   
 * */
#include <math.h> 
#include <stdio.h> 
#include <string.h> 
#include <stdarg.h>
#include <random>
#include "option.h"
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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


/* Generate a random number between 0 and 1 (excluded) from a uniform distribution. */
const gsl_rng_type * T;
gsl_rng * r;
double ran1(){
    double u = gsl_rng_uniform(r);
    // if (u > RNMX) ran1();
    return u;
}

/* Generate a random number from a Gaussian distribution of mean 0, and std 
   deviation 1.0. */
double gasdev(){
    return gsl_ran_ugaussian(r);
}

// --- define global parameters ------
static int ncomp = 10; // 7xthin + thick + bar, V, I, J, H, Ks 
static double tSFR = 7.0;  // time scale of SFR, 7.0 Gyr
static double rhot0;
static double MiniWDmax= 9; // 

// --- for fit to tE --- (from get_chi2_for_tE.c)
static int agesD[250], agesB[50], agesND[10];
static double MinidieD[250], MinidieB[50], MinidieND[10];
static int nageD=0, nageB=0, nageND=0;
static double mageB = 9, sageB = 1, mageND = 7, sageND = 1; // in Gyr

// --- for Mass function ---
static int nm;
static double logMst, dlogM;

// --- Parameters for bulge ---
// Values here will be overwritten in store_IMF_nBs using given IMF
static double fb_MS   = 1.62/2.07; // MS mass / total mass in bulge / bar
static double m2nb_MS  = 1/0.227943; // Msun/star in bulge / bar
static double m2nb_WD  = 1/0.847318; // Msun/WD   in bulge / bar
static double nMS2nRGb = 2.33232e-03; // n_RG/n_MS for bulge / bar
static double rho0b, n0MSb, n0RGb, n0b;

// Nuclear disk (for |b| < 1 deg.)
static int ND = 0, x0ND = 250, y0ND = 125, z0ND = 50;
static double C1ND = 2, rho0ND, n0MSND, n0RGND, n0ND;
static double fND_MS    = 0; // MS mass / total mass in NSD
static double m2nND_MS  = 0; // Msun/star in NSD
static double m2nND_WD  = 0; // Msun/WD   in NSD
static double nMS2nRGND = 0; // n_RG/n_MS for NSD

// --- Parameters for Disk ---
// Density values will be overwritten in store_IMF_nBs using given IMF
static double rho0d[8]   = {5.16e-03+3.10e-04, 5.00e-03+5.09e-04, 3.85e-03+5.42e-04, 3.18e-03+5.54e-04,
                            5.84e-03+1.21e-03, 6.24e-03+1.51e-03, 1.27e-02+3.49e-03, 1.68e-03+6.02e-04};
static double n0d[8]     = {1.51e-02+1.12e-04, 1.66e-02+3.22e-04, 1.40e-02+4.39e-04, 1.22e-02+5.15e-04, 
                            2.36e-02+1.25e-03, 2.63e-02+1.67e-03, 5.55e-02+4.08e-03, 7.91e-03+7.81e-04};
static double n0MSd[8]   = {1.51e-02, 1.66e-02, 1.40e-02, 1.22e-02, 2.36e-02, 2.63e-02, 5.55e-02, 7.91e-03};
static double n0RGd[8]   = {7.09e-06, 3.40e-05, 4.32e-05, 2.16e-05, 6.60e-05, 6.19e-05, 1.29e-04, 9.38e-06};
// Scale lengths and heights are fixed
static double y0d[3];
static int       Rd[3] = {5000, 2600, 2200};
static int        Rh = 3740, Rdbreak = 5300, nh = 1;
static double zd[8]   =  {61.47, 141.84, 224.26, 292.36, 372.85, 440.71, 445.37, 903.12};
static double zd45[8] =  {36.88,  85.10, 134.55, 175.41, 223.71, 264.42, 267.22, 903.12};
static int DISK, hDISK, addX, model;
static double R0, thetaD, x0_1, y0_1, z0_1=0, C1, C2, C3, Rc, frho0b, costheta, sintheta, zb_c;
static double x0_X, y0_X, z0_X=0, C1_X, C2_X, b_zX, fX, Rsin, b_zY, Rc_X;

//--- To give coordinate globally ---
static double *lDs, *bDs;

//--- For read Data from LFeachBD.dat & inputs/NbleNall_bin.dat ----
//--- For rough source mag and color constraint ----
static int nMIs, nVIs;
static double *MIs, **CumuN_MIs, dILF;
static double *VIs, ***f_VI_Is, dVILF;

//--- For Circular Velocity ------
static int nVcs=0;
static double Rcs[60], Vcs[60];

//--- Sun kinematics ------
static double vxsun = -10.0, Vsun = 11.0, vzsun = 7.0, vysun = 243.0;

//--- For Disk kinematics ------
static double ****fgsShu, ****PRRgShus, ****cumu_PRRgs;
static int ***n_fgsShu, ****kptiles;
static double hsigUt, hsigWt, hsigUT, hsigWT, betaU, betaW, sigU10d, sigW10d, sigU0td, sigW0td;
static double medtauds[8] = {0.075273, 0.586449, 1.516357, 2.516884, 4.068387, 6.069263, 8.656024, 12};
/* The line of sight toward (lSIMU, bSIMU) until Dmax pc needs to be inside the cylinder defined by
 * R < RenShu and -zenShu < z < zenShu 
 * Please change the following zenShu and/or RenShu value when you want to extend 
 * the line of sight outside of the default cylinder. */
static int zstShu =   0, zenShu = 3600, dzShu = 200;
static int RstShu = 500, RenShu = 9200, dRShu = 100; // use value @ RstShu for R < RstShu

//--- For Bulge kinematics ------
static int model_vb, model_vbz;
static double Omega_p, x0_vb, y0_vb, z0_vb, C1_vb, C2_vb, C3_vb, sigx_vb, sigy_vb, sigz_vb, vx_str, y0_str;
static double sigx_vb0, sigy_vb0, sigz_vb0;
static double x0_vbz, y0_vbz, z0_vbz, C1_vbz, C2_vbz, C3_vbz;

//--- For NSD, to store values of input_files/NSD_moments.dat ------
static double **logrhoNDs, **vphiNDs, ***logsigvNDs, **corRzNDs;
static double zstND = 0, zenND =  400, dzND = 5;
static double RstND = 0, RenND = 1000, dRND = 5;
static int nzND, nRND;

//--- Parameters to put Sgr A* on the GC ------
static double xyzSgrA[3] = {};

// Declare functions
double getx2y_khi(int n, double *x, double *y, double xin, int *khi);
double getx2y_ist(int n, double *x, double *y, double xin, int *ist);
double interp_x(int n, double *F, double xst, double dx, double xreq);
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq);
void   interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq);
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);

int main(int argc,char **argv)
{
  //--- read parameters ---
  int CheckD   = getOptiond(argc,argv,"CheckD", 1, 0);
  long seed    = getOptioni(argc,argv,"seed", 1, 12304357); // seed of random number
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, seed); 
  //--- Set params for Galactic model (default: E+E_X model in Koshimoto+2021) ---
  double M0_B      = getOptiond(argc,argv,"M0", 1, 1.0);
  double M1_B      = getOptiond(argc,argv,"M1", 1, 0.859770466578045);
  double M2_B      = getOptiond(argc,argv,"M2", 1, 0.08);
  double M3_B      = getOptiond(argc,argv,"M3", 1, 0.01);
  double Ml        = getOptiond(argc,argv,"Ml", 1, 0.001); // default : w/o planetary mass
  double Mu        = getOptiond(argc,argv,"Mu", 1, 120); // need to be fixed!!!. Affect normalizing bulge coeffs 
  double alpha1_B  = getOptiond(argc,argv,"alpha1", 1, -2.32279457078378);
  double alpha2_B  = getOptiond(argc,argv,"alpha2", 1, -1.13449983242887);
  double alpha3_B  = getOptiond(argc,argv,"alpha3", 1, -0.175862190587576);
  double alpha0_B  = getOptiond(argc,argv,"alpha0", 1,  alpha1_B);
  double alpha4_B  = getOptiond(argc,argv,"alpha4", 1,  alpha3_B);
  DISK     = getOptiond(argc,argv,"DISK",   1,    2); // 0: wo disk, 1: w/ disk+hole, 2: w/ disk like P17
  rhot0    = getOptiond(argc,argv,"rhot0",   1,   0.042); // local thin disk density, Msun/pc^3 (Bovy17: 0.042 +- 0.002 incl.BD)
  hDISK     = getOptiond(argc,argv,"hDISK",  1,    0); // 0: const scale height, 1: linear scale height
  addX      = getOptiond(argc,argv,"addX",   1,    5); // 0: no X-shape,  >=5: use model==addX as X-shape 
  model     = getOptiond(argc,argv,"model",  1,    5); // 
  R0     = getOptiond(argc,argv,"R0",     1,   8160); // 
  thetaD = getOptiond(argc,argv,"thetaD", 1,     27); // 
  frho0b = getOptiond(argc,argv,"frho0b", 1,  0.839014514507754); // 
  Rc     = getOptiond(argc,argv,"Rc", 1,  2631.78535429573); //
  zb_c   = getOptiond(argc,argv,"zb_c", 1,  1e+6); //
  if (model >= 4 && model <= 8){
    x0_1     = getOptiond(argc,argv,"x0", 1,  930.623146993329); // 
    y0_1     = getOptiond(argc,argv,"y0", 1,  370.784386649364); // 
    z0_1     = getOptiond(argc,argv,"z0", 1,  239.547516030578); // 
    C1     = getOptiond(argc,argv,"C1", 1,  1.20011972384328); // 
    C2     = getOptiond(argc,argv,"C2", 1,  4.09326795684828); // 
    C3     = getOptiond(argc,argv,"C3", 1,  1.0000); // 
  }
  if (addX >= 5){
    x0_X = getOptiond(argc,argv,"x0_X", 1,  278.027059842233); // 
    y0_X = getOptiond(argc,argv,"y0_X", 1,  176.318528789193); // 
    z0_X = getOptiond(argc,argv,"z0_X", 1,  286.791941602401); // 
    C1_X = getOptiond(argc,argv,"C1_X", 1,  1.3087131258784); // 
    C2_X = getOptiond(argc,argv,"C2_X", 1,  2.21745322869032); // 
    b_zX = getOptiond(argc,argv,"b_zX", 1,  1.37774815817195); // b_zX, slope of "X" of X-shape
    fX   = getOptiond(argc,argv,"fX",   1,  1.43975636704683); // fraction of X-shape
    Rc_X = getOptiond(argc,argv,"Rc_X",  1,  1301.63829617294); // 
  }
  b_zY   = getOptiond(argc,argv,"b_zY", 1, 0); //

  // ----- Kinematic parameters ------
  // for bar kinematic
  Omega_p  = getOptiond(argc,argv,"Omega_p",  1, 47.4105844018699);
  model_vb = getOptiond(argc,argv,"model_vb", 1,    5); // 
  x0_vb    = getOptiond(argc,argv,"x0_vb"  ,  1, 858.106595717275);
  y0_vb    = getOptiond(argc,argv,"y0_vb"  ,  1, 3217.04987721548);
  z0_vb    = getOptiond(argc,argv,"z0_vb"  ,  1, 950.690583433628);
  C1_vb    = getOptiond(argc,argv,"C1_vb"  ,  1, 4.25236641149869);
  C2_vb    = getOptiond(argc,argv,"C2_vb"  ,  1, 1.02531652066343);
  C3_vb    = getOptiond(argc,argv,"C3_vb"  ,  1, 1);
  sigx_vb  = getOptiond(argc,argv,"sigx_vb",  1, 151.854794853683);
  sigy_vb  = getOptiond(argc,argv,"sigy_vb",  1, 78.0278905748233);
  sigz_vb  = getOptiond(argc,argv,"sigz_vb",  1, 81.9641955092164);
  sigx_vb0 = getOptiond(argc,argv,"sigx_vb0",  1,  63.9939241108675);
  sigy_vb0 = getOptiond(argc,argv,"sigy_vb0",  1,  75.8180486866697);
  sigz_vb0 = getOptiond(argc,argv,"sigz_vb0",  1,  71.2336430487113);
  vx_str   = getOptiond(argc,argv,"vx_str" ,  1,    43.0364707040617);
  y0_str   = getOptiond(argc,argv,"y0_str" ,  1,    406.558313420815);
  model_vbz = getOptiond(argc,argv,"model_vbz",  1,    5); // 
  x0_vbz    = getOptiond(argc,argv,"x0_vbz"  ,  1, 558.430182718529);
  y0_vbz    = getOptiond(argc,argv,"y0_vbz"  ,  1, 2003.21703656302);
  z0_vbz    = getOptiond(argc,argv,"z0_vbz"  ,  1, 3823.20855045157);
  C1_vbz    = getOptiond(argc,argv,"C1_vbz"  ,  1, 3.71001266000693);
  C2_vbz    = getOptiond(argc,argv,"C2_vbz"  ,  1, 1.07455173734341);
  C3_vbz    = getOptiond(argc,argv,"C3_vbz"  ,  1, 1);

  // for disk kinematic (default: all-z + flat z_d^{thin} model in Koshimoto+21)
  hsigUt    = getOptiond(argc,argv,"hsigUt"  ,  1,  14300);// scale len of velo disp R (sigU) for thin
  hsigWt    = getOptiond(argc,argv,"hsigWt"  ,  1,   5900);// scale len of velo disp Z (sigW) for thin
  hsigUT    = getOptiond(argc,argv,"hsigUT"  ,  1, 180000);// scale len of velo disp R (sigU) for thick
  hsigWT    = getOptiond(argc,argv,"hsigWT"  ,  1, 9400);  // scale len of velo disp Z (sigW) for thick
  betaU     = getOptiond(argc,argv,"betaU"   ,  1, 0.32);  //  slope of age-sigU for thin
  betaW     = getOptiond(argc,argv,"betaW"   ,  1, 0.77);  //  slope of age-sigW for thin
  sigU10d   = getOptiond(argc,argv,"sigU10d" ,  1, 42.0);  // sigU for 10Gyr thin @Sunposi 
  sigW10d   = getOptiond(argc,argv,"sigW10d" ,  1, 24.4);  // sigW for 10Gyr thin @Sunposi
  sigU0td   = getOptiond(argc,argv,"sigU0td" ,  1, 75.0);  // sigU for thick @Sunposi
  sigW0td   = getOptiond(argc,argv,"sigW0td" ,  1, 49.2);  // sigW for thick @Sunposi

  // Use one of named models in Koshimoto+21
  int E_fg0 = getOptiond(argc,argv,"E_fg0", 1, 0);
  int G_fg0 = getOptiond(argc,argv,"G_fg0", 1, 0);
  int EXE_fg0 = getOptiond(argc,argv,"EXE_fg0", 1, 0);
  int GXG_fg0 = getOptiond(argc,argv,"GXG_fg0", 1, 0);
  if (E_fg0 == 1){ // E model
    model = 5, addX = 0;
    M0_B = 1.0, M1_B = 0.843651488650385, M2_B = 0.08, M3_B = 0.01;
    alpha1_B = -2.30708461042964, alpha2_B = -1.09811414023325, alpha3_B = -0.176687444667866;
    alpha0_B = alpha1_B, alpha4_B = alpha3_B;
    R0= 8160, thetaD = 27, 
    frho0b = 0.847695765083198, Rc = 2804.94024639663;
    x0_1 = 668.323640191308, y0_1 = 277.674592258175, z0_1 = 235.344943180979, 
    C1 = 1.40903573470129, C2 = 3.3497118832179, C3 = 1;
    model_vb = 5, model_vbz = 5;
    Omega_p = 49.5149910609312, vx_str = 48.7482280102778, y0_str = 392.515724264323,
    sigx_vb  = 156.055410564041, sigy_vb  = 83.8197043324931, sigz_vb  = 86.3564038759999,
    sigx_vb0 = 63.8292191277825, sigy_vb0 = 74.9469462226124, sigz_vb0 = 72.3085487545662,
    x0_vb  = 823.387929122523, y0_vb  = 9288.51482678556, z0_vb  = 864.479916419292,
    C1_vb  = 3.82820123451928, C2_vb  = 1.00573720627546,
    x0_vbz = 511.063328964278, y0_vbz = 2896.01606378595, z0_vbz = 2189.7664883434,
    C1_vbz = 3.04214421342047, C2_vbz = 1.00609904766722;
  }
  if (G_fg0 == 1){ // G model
    model = 6, addX = 0;
    M0_B = 1.0, M1_B = 0.896557393600988, M2_B = 0.08, M3_B = 0.01;
    alpha1_B = -2.39628188518525, alpha2_B = -1.18451896148506, alpha3_B = 0.168672130848533;
    alpha0_B = alpha1_B, alpha4_B = alpha3_B;
    R0= 8160, thetaD = 27, 
    frho0b = 0.777347874844233, Rc = 4838.85613149588;
    x0_1 = 1025.42128394916, y0_1 = 457.419718281149, z0_1 = 396.048253079423, 
    C1 = 2.00928445577057, C2 = 3.9678518191928, C3 = 1;
    model_vb = 5, model_vbz = 5;
    Omega_p = 40.5174879673548, vx_str = 11.9026090372449, y0_str = 20.1384817812277,
    sigx_vb  = 136.435675357212, sigy_vb  = 109.313291840218, sigz_vb  = 101.291432907346,
    sigx_vb0 = 76.0453005937702, sigy_vb0 = 67.9783132842431, sigz_vb0 = 74.7117386554542,
    x0_vb  = 1031.18302251324, y0_vb  = 2145.45565210108, z0_vb  = 727.233943973984,
    C1_vb  = 4.9302429910108, C2_vb  = 1.04038121792228,
    x0_vbz = 517.854475368706, y0_vbz = 1436.21008855387, z0_vbz = 1095.79181359292,
    C1_vbz = 2.3091601785779, C2_vbz = 1.03832670354301;
  }
  if (EXE_fg0 == 1){ // E+E_X model
    model = 5, addX = 5;
    M0_B = 1.0, M1_B = 0.859770466578045, M2_B = 0.08, M3_B = 0.01;
    alpha1_B = -2.32279457078378, alpha2_B = -1.13449983242887, alpha3_B = -0.175862190587576;
    alpha0_B = alpha1_B, alpha4_B = alpha3_B;
    R0= 8160, thetaD = 27, 
    frho0b = 0.839014514507754, Rc = 2631.78535429573;
    x0_1 = 930.623146993329, y0_1 = 370.784386649364, z0_1 = 239.547516030578, 
    C1 = 1.20011972384328, C2 = 4.09326795684828, C3 = 1;
    model_vb = 5, model_vbz = 5;
    Omega_p = 47.4105844018699, vx_str = 43.0364707040617, y0_str = 406.558313420815,
    sigx_vb  = 151.854794853683, sigy_vb  = 78.0278905748233, sigz_vb  = 81.9641955092164,
    sigx_vb0 = 63.9939241108675, sigy_vb0 = 75.8180486866697, sigz_vb0 = 71.2336430487113,
    x0_vb  = 858.106595717275, y0_vb  = 3217.04987721548, z0_vb  = 950.690583433628,
    C1_vb  = 4.25236641149869, C2_vb  = 1.02531652066343,
    x0_vbz = 558.430182718529, y0_vbz = 2003.21703656302, z0_vbz = 3823.20855045157,
    C1_vbz = 3.71001266000693, C2_vbz = 1.07455173734341;
    x0_X = 278.027059842233, y0_X = 176.318528789193, z0_X = 286.791941602401,
    C1_X = 1.3087131258784, C2_X = 2.21745322869032, 
    b_zX = 1.37774815817195, fX = 1.43975636704683, Rc_X = 1301.63829617294;
  }
  if (GXG_fg0 == 1){ // G+G_X model
    model = 6, addX = 6;
    M0_B = 1.0, M1_B = 0.901747918318042, M2_B = 0.08, M3_B = 0.01;
    alpha1_B = -2.32055781291126, alpha2_B = -1.16146692073597, alpha3_B = -0.222751835826612;
    alpha0_B = alpha1_B, alpha4_B = alpha3_B;
    R0= 8160, thetaD = 27, 
    frho0b = 0.861982105059042, Rc = 2834.43172768484;
    x0_1 = 1564.78976595399, y0_1 = 721.729645984158, z0_1 = 494.669973292979, 
    C1 = 1.20141097225, C2 = 3.09254667088709, C3 = 1;
    model_vb = 5, model_vbz = 5;
    Omega_p = 45.9061365175252, vx_str = 28.250608437116, y0_str = 11.4387290790323,
    sigx_vb  = 154.984185643613, sigy_vb  = 78.4783157632334, sigz_vb  = 83.2424209150283,
    sigx_vb0 = 63.3834790223473, sigy_vb0 = 75.1951371572303, sigz_vb0 = 69.6076680158332,
    x0_vb  = 939.470002303028, y0_vb  = 4228.61947632437, z0_vb  = 883.716365308057,
    C1_vb  = 4.59067123072475, C2_vb  = 1.00961963171066,
    x0_vbz = 699.073733500672, y0_vbz = 1729.91970395558, z0_vbz = 2028.24030134845,
    C1_vbz = 4.84589813971794, C2_vbz = 1.01718557457505;
    x0_X = 755.975821023038, y0_X = 312.17136920671, z0_X = 399.287597819655,
    C1_X = 1.21131134854495, C2_X = 1.30388556329566,
    b_zX = 1.37711800325276, fX = 2.99985800759016, Rc_X = 5174.00544959931;
  }

  costheta = cos(thetaD/180.0*PI) , sintheta = sin(thetaD/180.0*PI);

  // To put Sgr A* on the GC
  int CenSgrA = getOptioni(argc,argv, "CenSgrA", 1, 1);
  double lSgrA = -0.056; // Bland-Hawthorn & Gerhard 2016
  double bSgrA = -0.046;
  if (CenSgrA == 1){
    Dlb2xyz(R0, lSgrA, bSgrA, R0, xyzSgrA);
  }

  // Store Mass Function and calculate normalization factors for density distributions
  void store_IMF_nBs(int B, double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles, double M0, double M1, double M2, double M3, double Ml, double Mu, double alpha1, double alpha2, double alpha3, double alpha4, double alpha0);
  nm = 1000;
  double *logMass_B, *PlogM_cum_norm_B, *PlogM_B;
  int *imptiles_B;
  logMass_B        = (double*)calloc(nm+1, sizeof(double *));
  PlogM_B          = (double*)calloc(nm+1, sizeof(double *));
  PlogM_cum_norm_B = (double*)calloc(nm+1, sizeof(double *));
  imptiles_B       = (int*)calloc(22, sizeof(int *));
  store_IMF_nBs(1, logMass_B, PlogM_B, PlogM_cum_norm_B, imptiles_B, M0_B, M1_B, M2_B, M3_B, Ml, Mu, alpha1_B, alpha2_B, alpha3_B, alpha4_B, alpha0_B);

  // Make LF or VI vs I for the source
  double lSIMU   = getOptiond(argc,argv,"l",  1,  1.0);
  double bSIMU   = getOptiond(argc,argv,"b",  1, -3.9);
  double Isst    = getOptiond(argc,argv,"Isrange", 1, 14.0); // for Ds dist., default 14 < Is < 21, need AIrc
  double Isen    = getOptiond(argc,argv,"Isrange", 2, 21.0); // for Ds dist., default 14 < Is < 21, need AIrc
  double VIsst   = getOptiond(argc,argv,"VIsrange", 1, 0.0); // for Ds dist. w/ source col. const.
  double VIsen   = getOptiond(argc,argv,"VIsrange", 2, 0.0); // for Ds dist. w/ source col. const.
  double AIrc    = getOptiond(argc,argv,"AIrc", 1, 0); // AI for RC
  double EVIrc   = getOptiond(argc,argv,"EVIrc", 1, 0); // E(V-I) for RC
  double DMrc    = getOptiond(argc,argv,"DMrc", 1, 0); //mean DM for RC, default is given later
  if (fabs(lSIMU) > 10 || fabs(bSIMU) > 7 || fabs(bSIMU) < 1.5){
    printf("# WARNING: genulens is designed for use in |l| < ~10 and ~1.5 < |b| < ~7 and (l,b)= (%.3f, %.3f) deg. is outside of the range.\n",lSIMU,bSIMU);
  }
  // When only Isrange is given
  int narry = 960;
  int make_LFs(double *MIs, double **CumuN_MIs, double *logMass, double *PlogM_cum_norm);
  if (Isen - Isst > 0 && VIsen - VIsst == 0 && AIrc > 0){
    CumuN_MIs = (double**)malloc(sizeof(double *) * ncomp);
    for (int i=0; i<ncomp; i++){
       CumuN_MIs[i] = (double*)calloc(narry, sizeof(double *));
    }
    MIs = (double*)malloc(sizeof(double *) * narry);
    nMIs = make_LFs(MIs, CumuN_MIs, logMass_B, PlogM_cum_norm_B);
    dILF = (MIs[nMIs-1] - MIs[0])/(nMIs - 1);
  }

  // When VIsrange is also given
  void store_VI_MI(double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI, double *MIs, double *VIs, double ***f_VI_Is, double *logMass, double *PlogM_cum_norm);
  if (Isen - Isst > 0 && VIsen - VIsst > 0 && AIrc > 0 && EVIrc > 0){
    double MIst = -5, MIen = 10, VIst = 0.0, VIen = 3.0;
    nMIs = 150, nVIs = 30; // dI = 0.1, dVI = 0.1
    MIs = (double*)calloc(nMIs+1, sizeof(double *));
    VIs = (double*)calloc(nVIs+1, sizeof(double *));
    f_VI_Is = (double***)malloc(sizeof(double *) * ncomp);
    for (int i=0; i<ncomp; i++){
      f_VI_Is[i] = (double**)malloc(sizeof(double *) * (nVIs+1));
      for (int j=0; j<nVIs+1; j++){
        f_VI_Is[i][j] = (double*)calloc(nMIs+1, sizeof(double *));
      }
    }
    store_VI_MI(MIst, MIen, nMIs, VIst, VIen, nVIs, MIs, VIs, f_VI_Is, logMass_B, PlogM_cum_norm_B);
    dILF  = (double) (MIen - MIst)/nMIs;
    dVILF = (double) (VIen - VIst)/nVIs;
  }

  // Read empirical mass-luminocity relation to calculate lens brightness
  int ncols = 5; // V_j, I_c, J_2MASS, H_2MASS, K_2MASS
  double *M_emps, **Mag_emps;
  narry = 60;
  M_emps = (double*)calloc(narry, sizeof(double *));
  Mag_emps = (double**)malloc(sizeof(double *) * ncols);
  for (int i=0; i<ncols; i++){
     Mag_emps[i] = (double*)calloc(narry, sizeof(double *));
  }
  char *file_MLemp = (char*)"input_files/MLemp.dat";
  int read_MLemp(char *infile, double *M_emps, double **Mag_emps);
  int nMLemp = read_MLemp(file_MLemp, M_emps, Mag_emps);

  // Store Cumu P_Shu
  int nfg = 100;
  int nz = (zenShu - zstShu)/dzShu + 1;
  int nR = (RenShu - RstShu)/dRShu + 1;
  int ndisk = 8;
  fgsShu      = (double****)malloc(sizeof(double *) * nz);
  PRRgShus    = (double****)malloc(sizeof(double *) * nz);
  cumu_PRRgs  = (double****)malloc(sizeof(double *) * nz);
  n_fgsShu  = (int***)malloc(sizeof(int *) * nz);
  kptiles   = (int****)malloc(sizeof(int *) * nz);
  for (int i=0; i<nz; i++){
    fgsShu[i]  = (double***)malloc(sizeof(double *) * nR);
    PRRgShus[i] = (double***)malloc(sizeof(double *) * nR);
    cumu_PRRgs[i] = (double***)malloc(sizeof(double *) * nR);
    n_fgsShu[i]  = (int**)malloc(sizeof(int *) * nR);
    kptiles[i]   = (int***)malloc(sizeof(int *) * nR);
    for (int j=0; j<nR; j++){
      fgsShu[i][j]  = (double**)malloc(sizeof(double *) * ndisk);
      PRRgShus[i][j] = (double**)malloc(sizeof(double *) * ndisk);
      cumu_PRRgs[i][j] = (double**)malloc(sizeof(double *) * ndisk);
      kptiles[i][j] = (int**)malloc(sizeof(int *) * ndisk);
      n_fgsShu[i][j]  = (int*)calloc(ndisk, sizeof(int *));
      for (int k=0; k<ndisk; k++){
        fgsShu[i][j][k]   = (double*)calloc(nfg, sizeof(double *));
        PRRgShus[i][j][k] = (double*)calloc(nfg, sizeof(double *));
        cumu_PRRgs[i][j][k] = (double*)calloc(nfg, sizeof(double *));
        kptiles[i][j][k]  = (int*)calloc(22, sizeof(int *));
      }
    }
  }
  char *fileVc = (char*)"input_files/Rotcurve_BG16.dat";
  void store_cumuP_Shu(char *infile);
  store_cumuP_Shu(fileVc);

  // set y0d for disk normalize
  y0d[0] = (DISK == 1) ? exp(-R0/Rd[0] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[0]);
  y0d[1] = (DISK == 1) ? exp(-R0/Rd[1] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[1]);
  y0d[2] = (DISK == 1) ? exp(-R0/Rd[2] - pow(((double)Rh/R0),nh))  :  exp(-R0/Rd[2]);

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
  printf("#            tau   Rd  zd zd45 sigU0 sigW0  RsigU  RsigW    rho0        n0     n0WD \n");
  printf("#            Gyr   pc  pc   pc  km/s  km/s     pc     pc  Msun/pc^3  */pc^3   */pc^3\n");
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
    printf ("#   Disk%d: %5.2f %4d %3.0f  %3d %5.2f %5.2f %6.0f %6.0f   %.2e %.2e %.2e\n",i+1,medtauds[i], rd, zd[i],zdtmp,sigU0,sigW0,hsigU,hsigW,rho0d[i],n0d[i],n0d[i]-n0MSd[i]);
  }

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
  
  // Calc PA (added on 2022/6/14)
  double PA, cosPA, sinPA;
  void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
  calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);
  
  // Read Input parameters for simulation
  long   NSIMU     = getOptiond(argc,argv,"Nsimu",  1, 100000);
  long   NlikeMIN  = getOptiond(argc,argv,"NlikeMIN",  1, 0);
  // double PA         = getOptiond(argc,argv,"PA", 1, 59.56); // postition angle l to E value. 59.56 is for (l,b) = (0.94, -1.48), needed to calculate piEn and piEe
  // double cosPA    = cos(PA/180.0*PI), sinPA = sin(PA/180.0*PI);
  // printf("PA= %f, cosPA= %f, sinPA= %f\n",PA,cosPA,sinPA);
  double vEarthl  = getOptiond(argc,argv,"vEarthlb", 1, 11.9392); // in km/s, default for MB16227 occured on May
  double vEarthb  = getOptiond(argc,argv,"vEarthlb", 2,-17.9020); // in km/s, default for MB16227 occured on May
  double vEarthE  = getOptiond(argc,argv,"vEarthEN", 1, 0); // in km/s
  double vEarthN  = getOptiond(argc,argv,"vEarthEN", 2, 0); // in km/s
  if (vEarthl == 11.9392 && vEarthb == -17.9020 && vEarthN != 0 && vEarthE != 0){
    vEarthb =  cosPA*vEarthN - sinPA*vEarthE;
    vEarthl =  sinPA*vEarthN + cosPA*vEarthE;
  }
  double gammaDs  = getOptiond(argc,argv,"gammaDs",  1, 0.5);
  double wtD_L    = getOptiond(argc,argv,"wtD_L",    1, 0);  // parameter for importance sampling 
  double wtM_L    = getOptiond(argc,argv,"wtM_L",    1, 0);  // parameter for importance sampling
  int SMALLGAMMA  = getOptiond(argc,argv,"SMALLGAMMA",   1,  0);
  int VERBOSITY   = getOptiond(argc,argv,"VERBOSITY",  1, 0);
  int UNIFORM     = getOptiond(argc,argv,"UNIFORM",   1,  0);
  int BINARY      = getOptiond(argc,argv,"BINARY",   1,  0);
  int REMNANT     = getOptiond(argc,argv,"REMNANT",  1,  0);
  int onlyWD      = getOptiond(argc,argv,"onlyWD",  1,  0);
  if (onlyWD == 1) REMNANT = 0;
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
  // double femusE  = getOptiond(argc,argv,"musE", 3, 0); // 
  double femusE = 0;
  int    musRCG  = getOptiond(argc,argv,"musRCG", 1, 0); // 0: mus relative to Sun, 1: mus relative to RCG 
  double ILobs = getOptiond(argc,argv,"IL", 1, 14.00); // mag
  double ILe   = getOptiond(argc,argv,"IL", 2,  0.01); // mag
  // double feIL  = getOptiond(argc,argv,"IL", 3, 0); //
  double feIL = 0; 
  int    ILdet = getOptiond(argc,argv,"ILdet", 1, 2); // 0: det, 1: upper limit, 2: lower limit, default: 2
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
  double AI0  = (hscale > 0 && Dmean > 0) ? AIrc  / (1 - exp(-Dmean/hscale)) : 0; // 
  double EVI0 = (hscale > 0 && Dmean > 0) ? EVIrc / (1 - exp(-Dmean/hscale)) : 0; // 
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
  printf("#    CenSgrA= %d     (0: GC at (l,b)=(0,0), 1: GC at (l,b)= (%.3f, %.3f))\n", CenSgrA, lSgrA, bSgrA);
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
  if (ILe > 0)    printf("#      IL = %.3f +- %.3f Fe_add= %6.1f det= %d (0: det, 1: upper limit, 2: lower limit)\n",ILobs,ILe,feIL,ILdet);
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
  double rhos[10] = {}, xyz[3] = {}, xyb[2] = {};
  // Lens   : include REMNANT, mass basis 
  // Source : only stars, number basis 
  int Dmax = 16000;
  // int Dmax = 12000;
  int nbin = (ND > 0 && fabs(lSIMU) < 0.05 && fabs(bSIMU) < 0.05) ? 0.20*Dmax+0.5 
           : (ND > 0 && fabs(lSIMU) < 0.10 && fabs(bSIMU) < 0.10) ? 0.10*Dmax+0.5
           : (ND > 0) ? 0.04*Dmax+0.5 
           : 0.01*Dmax+0.5;
  double dD = (double) Dmax/nbin;
  double *D, **rhoD_S, **rhoD_L, **cumu_rho_S, **cumu_rho_L, *cumu_rho_all_S, *cumu_rho_all_L;
  D               = (double *)calloc(nbin+1, sizeof(double *));
  cumu_rho_all_S  = (double *)calloc(nbin+1, sizeof(double *));
  cumu_rho_all_L  = (double *)calloc(nbin+1, sizeof(double *));
  rhoD_S      = (double **)malloc(sizeof(double *) * ncomp);
  rhoD_L      = (double **)malloc(sizeof(double *) * ncomp);
  cumu_rho_S  = (double **)malloc(sizeof(double *) * ncomp);
  cumu_rho_L  = (double **)malloc(sizeof(double *) * ncomp);
  for (int i=0; i<ncomp; i++){
    rhoD_S[i]     = (double *)calloc(nbin+1, sizeof(double *));
    rhoD_L[i]     = (double *)calloc(nbin+1, sizeof(double *));
    cumu_rho_S[i] = (double *)calloc(nbin+1, sizeof(double *));
    cumu_rho_L[i] = (double *)calloc(nbin+1, sizeof(double *));
  }
  double fLF_detect(double extI, double Imin, double Imax, int idisk);
  double fIVI_detect(double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk);
  printf("#----- Mass density distribution along (l, b)= (%.3f, %.3f) w/ wtD_L= %.1f --------\n",lSIMU,bSIMU,wtD_L);
  int npri = (ND > 0) ? 40 : 10;
  double nallS = 0;
  for (int ibin=0; ibin<=nbin; ibin++){
    D[ibin] = (double) ibin/nbin * Dmax;
    calc_rho_each(D[ibin], idata, rhos, xyz, xyb);
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    if (ibin%npri==0) printf ("# %5.0f %5.0f %5.0f ",D[ibin],R,xyz[2]);
    double rhosum = 0;
    double extI  =  AI0 * (1 - exp(-D[ibin]/hscale)) + 5 * log10(0.1*(D[ibin] + 0.1));
    double extVI = EVI0 * (1 - exp(-D[ibin]/hscale));
    // printf ("%5.0f %7.3f ",D[ibin],extI);
    for (int i=0;i<ncomp;i++){
      double nMS = (i == 8) ? n0MSb*rhos[8] : (i == 9) ? n0MSND*rhos[9] : n0MSd[i]*rhos[i];
      double rho = (i == 8) ? n0b  *rhos[8] : (i == 9) ? n0ND  *rhos[9] : n0d[i]  *rhos[i];
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
      rhoD_L[i][ibin] = rho;
      // Cumulative. The min and max edges are not correctly treated. 
      // But practically OK especially considering a risk of ran < Cumu_min or ran > Cumu_max when strictly considered.
      cumu_rho_S[i][ibin]  = (ibin==0) ? 0 : cumu_rho_S[i][ibin-1] + 0.5*(rhoD_S[i][ibin-1] + rhoD_S[i][ibin])*dD;
      cumu_rho_L[i][ibin]  = (ibin==0) ? 0 : cumu_rho_L[i][ibin-1] + 0.5*(rhoD_L[i][ibin-1] + rhoD_L[i][ibin])*dD;
      cumu_rho_all_S[ibin] += cumu_rho_S[i][ibin];
      cumu_rho_all_L[ibin] += cumu_rho_L[i][ibin];
      rhosum += rho;
      if (ibin%npri==0){ 
        printf (" %d: %.1e ",i,rhoD_L[i][ibin]);
        printf ("( %.2e )",cumu_rho_L[i][ibin]);
      }
    }
    // printf ("\n");
    if (ibin%npri==0){ 
        printf (" All: %.1e ",rhosum);
        printf ("( %.2e )\n",cumu_rho_all_L[ibin]);
    }
  }
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
      printf ("%d %6.0f %d %6.0f\n",j_S,d_S,j_L,d_L);
      // pick velocities
      void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD); //
      double vxyz_L[3] = {};
      double tau_l = (j_L == 9) ? mageND + sageND*gasdev()
                   : (j_L == 8) ? mageB + sageB*gasdev() 
                   : medtauds[j_L];
      get_vxyz_ran(vxyz_L, j_L, tau_l, d_L, lDs[idata], bDs[idata]);
      double vx_l = vxyz_L[0];
      double vy_l = vxyz_L[1];
      double vz_l = vxyz_L[2];

      printf ("%d %6.0f %d %6.0f %7.2f %7.2f %7.2f\n",j_S,d_S,j_L,d_L, vx_l, vy_l, vz_l);
    }
    exit(1);
  }

  // release luminosity functions
  if (Isen - Isst > 0 && VIsen - VIsst == 0 && AIrc > 0){
    for (int i=0; i<ncomp; i++){
       free(CumuN_MIs[i]);
    }
    free(CumuN_MIs);
    free(MIs);
  }
  if (Isen - Isst > 0 && VIsen - VIsst > 0 && AIrc > 0 && EVIrc > 0){
    for (int i=0; i<ncomp; i++){
      for (int j=0; j<nVIs+1; j++){
        free(f_VI_Is[i][j]);
      }
      free(f_VI_Is[i]);
    }
    free(f_VI_Is);
    free(MIs);
    free(VIs);
  }

  // Monte Carlo simulation
  printf ("#----- Output of Monte Carlo simulation w/ VERBOSITY= %d and seed= %ld -------- \n",VERBOSITY,seed);
  double getcumu2xist (int n, double *x, double *F, double *f, double Freq, int ist, int inv);
  double ncntall = 0, ncnts = 0, ncntbWD = 0, ncntbCD = 0; 
  double nBD = 0, nMS = 0, nWD = 0, nNS= 0, nBH =  0;
  int Nlike = 0;
  long NrejIS = 0, Ngen = 0;
  double wtlike = 0, wtlike_tE =0;
  double SumGamma = 0, SumtE = 0; // for mean tE
  double logtEmin = -1, logtEmax = 2, NbintE = 300, NlogtEs[500] = {}; // for median tE
  if (VERBOSITY == 2) printf ("#        wtj           tE       thetaE          piEN          piEE   D_S         muSl         muSb iS iL fREM");
  if (VERBOSITY == 3) printf ("#        wtj          M_L   D_L   D_S          t_E      theta_E         pi_E         pi_EN         pi_EE       mu_rel        mu_Sl        mu_Sb     I_L iS iL fREM");
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
         double vkick = (iREM == 1) ? 350 : 100;  // Table 2 of Lam et al. 2020
         vxadd =  vkick * cos(thetakick) * cos(phikick);
         vyadd =  vkick * cos(thetakick) * sin(phikick);
         vzadd =  vkick * sin(thetakick);
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
     // double Gamma = 8e-09 * D_l*D_l*thetaE*murel; // 8e-09 makes Gamma to be < ~1
     // printf("Ml= %.6f, D_l= %.2f Gamma= %.5e piE= %.5f thetaE= %.5f tE= %.3f\n",M_l, D_l,Gamma,piE,thetaE,tE);
     // for mean tE
     SumGamma += Gamma*addGammaIS;
     SumtE    += Gamma*addGammaIS*tE;
     // for median tE
     int ilogtE = (log10(tE) - logtEmin)/(logtEmax - logtEmin)*NbintE;
     if (ilogtE < 0) ilogtE = 0;
     if (ilogtE > NbintE - 1) ilogtE = NbintE - 1;
     NlogtEs[ilogtE] += Gamma*addGammaIS;
     if (Gamma < ran1() && SMALLGAMMA == 0){
       j--;
       continue;
     }
     double like_obs(double mod, double obs, double err, double fe, int det, int UNIFORM);
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

     // Constraint from lens brigtness (currently only I)
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
       IL += 5 * log10(0.1*D_l) + AI0 * (1 - exp(-D_l/hscale)); // add DM and AI
       if(ILe > 0){  // This is true when default to cut too bright lens
         double Gamma_IL = like_obs(IL, ILobs, ILe, feIL, ILdet, UNIFORM);
         like_IL = (Gamma_IL > 0) ? 1 : 0;
         if (Gamma_IL > 0) addGamma *= Gamma_IL;
       }
     }

     // Combine all constraints
     Gamma *= addGamma * addGammaIS;
     wtj   *= addGamma * addGammaIS;
     int like = like_tE*like_thetaE*like_piE*like_mus*like_IL;
     // if (VERBOSITY == 1) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %2d %2d %d %d\n", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, i_s,i_l,fREM,like);
     if (like_tE == 1) wtlike_tE += wtj;
     if (like    == 1) Nlike ++, wtlike    += wtj;
     else continue;
     // if (VERBOSITY == 1) printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %1d %1d %5.2f %5.2f %d %.5e", wtj,tE,thetaE,piE,M_l,D_s, D_l,vx_s,vy_s,vz_s,vx_l,vy_l,vz_l,murel,i_s,i_l,tau_s,tau_l,fREM,Gamma);
     if (VERBOSITY == 2) printf("%12.6e %12.6e %12.6e %13.6e %13.6e %5.0f %12.5e %12.5e %2d %2d %d", wtj, tE,thetaE,piEN,piEE,D_s,muSl,muSb,i_s,i_l,fREM);
     if (VERBOSITY == 3) printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e %12.6e %12.5e %12.5e %7.3f %2d %2d %d", wtj, M_l, D_l, D_s, tE, thetaE, piE, piEN, piEE, murel, muSl,muSb, IL, i_s,i_l,fREM);
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
  gsl_rng_free(r);
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
  free(logMass_B       );
  free(PlogM_cum_norm_B);
  free(PlogM_B         );
  free(imptiles_B      );
  free(lDs);
  free(bDs);
  for (int i=0; i<nz; i++){
    for (int j=0; j<nR; j++){
      for (int k=0; k<ndisk; k++){
        free(fgsShu[i][j][k]);
        free(PRRgShus[i][j][k]);
        free(cumu_PRRgs[i][j][k]);
        free(kptiles[i][j][k]);
      }
      free(fgsShu[i][j]);
      free(PRRgShus[i][j]);
      free(cumu_PRRgs[i][j]);
      free(kptiles[i][j]);
      free(n_fgsShu[i][j]);
    }
    free(fgsShu[i]);
    free(PRRgShus[i]);
    free(cumu_PRRgs[i]);
    free(kptiles[i]);
    free(n_fgsShu[i]);
  }
  free(fgsShu);
  free(PRRgShus);
  free(cumu_PRRgs);
  free(kptiles);
  free(n_fgsShu);
  return 0;
} // end main

//----------------
void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA){
  // Calculate position angle. Added on 2022/6/14 copied from mu_Hlb2mu_Gne.pl
  double lNP = 122.9320; // galactic l of equatorial north pole, http://ned.ipac.caltech.edu/forms/calculator.html
  double bNP =  27.1284; // galactic b of equatorial north pole, http://ned.ipac.caltech.edu/forms/calculator.html
  double cosl   = cos(gl/180. * PI);
  double sinl   = sin(gl/180. * PI);
  double cosb   = cos(gb/180. * PI);
  double sinb   = sin(gb/180. * PI);
  double coslNP = cos(lNP/180. * PI);
  double sinlNP = sin(lNP/180. * PI);
  double cosbNP = cos(bNP/180. * PI);
  double sinbNP = sin(bNP/180. * PI);
  double nvector[3] = {cosb*cosl,cosb*sinl,sinb}; // unit vector along line of sight 
  double galNP[3]   = {0,0,1};                    // unit vector toward galactic north pole
  double eqNP[3] = {cosbNP*coslNP,cosbNP*sinlNP,sinbNP}; // unit vector toward equatorial north pole

  // Make Unit Vectors along (b, l) and along (N, E) on the target sky
  double dot(double *a, double *b);
  void cross(double *c, double *a, double *b);
  void norm_vec(double *a);
  double elvector[3] = {}, ebvector[3] = {}, eEvector[3] = {}, eNvector[3] = {};
  cross(elvector, galNP, nvector); // p x n
  norm_vec(elvector); // normalize
  cross(ebvector, nvector, elvector); // eb = n x el

  cross(eEvector, eqNP, nvector); // p x n
  norm_vec(eEvector); // normalize
  cross(eNvector, nvector, eEvector); // eN = n x eE

  // Calc -PA (negative position angle, from l to E or from b to N. A positive PA should be defined from N to b, maybe
  double crosstmp[3] = {};
  *cosPA = dot(elvector,eEvector);
  cross(crosstmp, elvector, eEvector);
  *sinPA = - dot(nvector, crosstmp);
  *PA = 180*atan2(*sinPA, *cosPA)/PI; // atan2(Y,X) -PA, radian, l to E
}

//----------------
void calc_opticaldepth(double *tauall, double *Nsall, int idata, int Dsmax21, double AI0, double hscale, double Isst, double Isen)
{
  /* For optical depth or event rate calculation    
   * Not optimized for this code and many calculations are dupulicated with previous calculations. */
  int get_p_integral(int nji, double *ls, double *ks);
  int nmin;
  // param for integration over Ds
  int Dsmin21   = 400;
  double dDs0 = 50;
  int nbun21  = (Dsmax21 - Dsmin21)/dDs0 + 0.5;
  int nji21   = 2;
  int narry21 = (nji21 <= 1) ?  1 :  (nji21 <= 2) ?  3 :  (nji21 <= 4) ?  9 :  
                (nji21 <= 6) ? 18 :  (nji21 <= 8) ? 30 : 42;
  double *ls21, *ks21;
  ls21 = (double *)malloc(sizeof(double *) * narry21);
  ks21 = (double *)malloc(sizeof(double *) * narry21);
  nmin = get_p_integral(nji21, ls21, ks21);
  if (nbun21 < nmin){
     printf ("# Warning: nbun21 (= %d) is updated to be nmin = %d\n",nbun21,nmin);
     nbun21 = nmin;
  }
  int ncalc21 = nbun21 + 1 + 2*narry21 - 2*nji21;  // needed number of array to conduct nji21 sekibun

  // param for integration over Dl
  int nbunDl0   = 20;
  int Dlmin     =  0;
  int njiDl     =  2;

  nbunDl0 = (nji21 == 2) ? (Dsmin21 - Dlmin)/(dDs0/2)
          : (nji21 == 1) ? (Dsmin21 - Dlmin)/dDs0 : nbunDl0; 
  double dDl0 = (double) (Dsmin21 - Dlmin  )/nbunDl0; // default:  4000/20 = 200
  int narryDl = (njiDl <= 1) ?  1 :  (njiDl <= 2) ?  3 :  (njiDl <= 4) ?  9 :  
                (njiDl <= 6) ? 18 :  (njiDl <= 8) ? 30 : 42;
  double *lsDl, *ksDl;
  lsDl = (double *)malloc(sizeof(double *) * narryDl);
  ksDl = (double *)malloc(sizeof(double *) * narryDl);
  nmin = get_p_integral(njiDl, lsDl, ksDl);
  if (nbunDl0 < nmin){ // nbunDl0 has to NOT be updated here
     printf ("Error: nbunDl0 (= %d) has to be larger than nmin (= %d)!!\n",nbunDl0,nmin);
     printf ("       Change njiDl or nbunDl0 so that this does not occur!!\n");
     exit(1);
  }
  int ncalcDl0 = nbunDl0 + 1 + 2*narryDl - 2*njiDl;  // needed number of array to conduct njiDl sekibun


  // Calc tau(Ds) for Dsmin21 <= Ds <= Dsmax21  (from calc_tauDs_nji2 in get_chi2_forN13_M19_C19_Neve_tE.c)
  // Prepare tauDs and NDs
  double *tauDs, *NDs, *wtDBs;
  tauDs = (double *)malloc(sizeof(double *) * ncalc21);  // tau[Ds]
  NDs   = (double *)malloc(sizeof(double *) * (ncalc21 + ncalcDl0)); // NDs stores ncalc21*NDs[Ds] + ncalcDl*NDs[Dl]
  wtDBs = (double *)malloc(sizeof(double *) * (ncalc21 + ncalcDl0)); // wtDBs:     0.1*int(wtDBs) -> mean diski
                                                           //      : wtDBs - int(wtDBs) -> disk/total 
  double calc_rho_n(double D, int idata, double *rho_n); // return rho(D) & n(D)
  double *rhoDlkpt0, *rhoDlkpt1;
  rhoDlkpt0 = (double *)calloc(narryDl, sizeof(double *));
  int narrytmp = (Dsmax21 - Dlmin)/dDl0; // default should be 16000/200 = 80 
  rhoDlkpt1 = (double *)calloc(narrytmp + 1, sizeof(double *));
  for (int iDs = 0; iDs < ncalc21; iDs++){    // loop for the next 21 integration
    int iDstmp = iDs - 2*narry21 + nji21;        // iDstmp = nji21 - (nbun21 - nji21)
    double Ds = (iDs>=2*narry21) ? Dsmin21 + dDs0 * iDstmp 
              : (iDs   % 2 == 0) ? Dsmin21 + dDs0 * ls21[iDs/2] : Dsmax21 - dDs0 * ls21[iDs/2];
    double Dlmax = Ds;
    int nbun = nbunDl0 + (Dlmax - Dsmin21)/dDl0; 
    double nbuntmp =   (Dlmax - Dlmin)/dDl0;
    if ((double) nbun != nbuntmp){ 
      printf("ERROR: nji21 (%d), dDl0 (%f), dDs0 (%f) has to satisfy the followings:\n",nji21,dDl0,dDs0);
      printf("       nji21 <= 2\n");
      printf("       dDl0 =     dDs0/n when nji21 == 1 (n: natural number)\n");
      printf("       dDl0 = 0.5*dDs0/n when nji21 == 2 (n: natural number)\n");
      exit(1);
    }
    tauDs[iDs] = 0;
    double rho_n[2] = {}, wtDB, Dl, tau, tau0; 
    for(int j=0;j< narryDl;j++){
        double dDltmp = dDl0*lsDl[j];
        // 以下はtau(Dlmin + dDl0*lsDl[j]) が計算されて、 Dlmin は iDs によらずに = 0 で、
        // 密度も同じなので、最初にrho*Dl にkeepしてそれを使いまわす
        Dl   =  Dlmin + dDltmp;
        if (iDs == 0){
          wtDB = calc_rho_n(Dl, idata, rho_n);
            NDs[ncalc21+2*j] = rho_n[1]*Dl*Dl; // store n[Ds], Ds = Dlmin - Dsmin21
          wtDBs[ncalc21+2*j] = wtDB; // store n[Ds], Ds = Dlmin - Dsmin21
          rhoDlkpt0[j] = rho_n[0] * Dl;
        }
        tau0 = rhoDlkpt0[j] * (1 - Dl/Ds); // tau = rho * Dl * (1 - Dl/Ds)

        // 以下はtau(Dlmax - dDl0*lsDl[j]) が計算されて、 Dlmax は更新されていくので、rho*Dlを使い回せない。
        // ただ、dDl0 の倍数のものに関しては後で違うiDsの時に使いまわせる可能性があるのでkeepする。
        Dl   =  Dlmax - dDltmp;
        if (ceil(lsDl[j]) == floor(lsDl[j]) && rhoDlkpt1[nbun-(int)lsDl[j]] > 0 && j != 0){
          tau = rhoDlkpt1[nbun-(int)lsDl[j]] * (1 - Dl/Ds);
          // printf (" iDs= %2d j= %2d lsDl= %.2f rhoDlkpt1[%d]= %.5e\n",iDs,j,lsDl[j],nbun-(int)lsDl[j],rhoDlkpt1[nbun-(int)lsDl[j]]);
        }else{
          wtDB = calc_rho_n(Dl, idata, rho_n);
          tau  = rho_n[0] * Dl * (1 - Dl/Ds);
          if (ceil(lsDl[j]) == floor(lsDl[j])) rhoDlkpt1[nbun-(int)lsDl[j]] = rho_n[0] * Dl;
          if (j   == 0)   NDs[iDs]         = rho_n[1]*Dl*Dl; // store n[Ds], Dsmin21 <= Ds <= Dsmax21
          if (j   == 0) wtDBs[iDs]         = wtDB;
          if (iDs == 0)   NDs[ncalc21+2*j+1] = rho_n[1]*Dl*Dl; // store n[Ds], Dlmin   <= Ds <= Dsmin21
          if (iDs == 0) wtDBs[ncalc21+2*j+1] = wtDB;
          // printf ("j= %2d lsDl= %.2f sig21[%2d]= %8.2f wtDB= %7.4f\n",j,lsDl[j],iDs,NDs[iDs]*dDs0*STR2MIN2,wtDBs[iDs]);
        }
        // 
        tauDs[iDs] += (tau0 + tau)*ksDl[j];
    }
    for(int j=njiDl;j<=nbun-njiDl;j++){
        Dl  =  Dlmin + dDl0*j;
        if (rhoDlkpt1[j] == 0){
          wtDB = calc_rho_n(Dl, idata, rho_n);
          if (iDs == 0)   NDs[ncalc21+2*narryDl+j-njiDl] = rho_n[1]*Dl*Dl; // store n[Ds], Dlmin <= Ds <= Dsmin21
          if (iDs == 0) wtDBs[ncalc21+2*narryDl+j-njiDl] = wtDB;
          rhoDlkpt1[j] = rho_n[0] * Dl;
        }
        tau = rhoDlkpt1[j] * (1 - Dl/Ds); // tau = rho * Dl * (1 - Dl/Ds)
        tauDs[iDs] += tau;
    }
    tauDs[iDs] *= dDl0;
    // printf ("%5d %8.4f %8.4f Ds[%02d]= %5.0f nbun= %d (dDl0= %6.2f) ls= %7.3f tauDs= %.4e sig21= %8.2f wtDB= %7.4f\n",idata,lDs[idata],bDs[idata],iDs,Ds,nbun,dDl0,(Ds-Dsmin21)/dDs0,tauDs[iDs]*PI4GC2,NDs[iDs]*dDs0*STR2MIN2,wtDBs[iDs]);
  }
  free(rhoDlkpt0);
  free(rhoDlkpt1);

  // Calc N_source and tau (from int_Ds21 in get_chi2_forN13_M19_C19_Neve_tE.c)
  double fLF_detect(double extI, double Imin, double Imax, int idisk);
  *tauall = 0, *Nsall = 0; 
  for (int iDs = 0; iDs < ncalc21; iDs++){
     int iDstmp = iDs - 2*narry21 + nji21;        // iDstmp = nji21 - (nbun21 - nji21)
     double Ds = (iDs>=2*narry21) ? Dsmin21 + dDs0 * iDstmp 
               : (iDs   % 2 == 0) ? Dsmin21 + dDs0 * ls21[iDs/2] : Dsmax21 - dDs0 * ls21[iDs/2];

     // Extinction + Distance modulus at Ds
     double extI =  AI0 * (1 - exp(-Ds/hscale)) + 5 * log10(0.1*(Ds + 0.1));

     // DISK fraction of 14 < I < 21
     double m_idisk = (int) wtDBs[iDs];
     double f_disk  = 10*(wtDBs[iDs] - m_idisk);
     m_idisk *= 0.1;
     double wtD = m_idisk - (int) m_idisk; // e.g.) 4.2 - 4.0 = 0.2 
     double fMagD =  (1 - wtD) * fLF_detect(extI, Isst, Isen, (int)floor(m_idisk))
                   +      wtD  * fLF_detect(extI, Isst, Isen, (int) ceil(m_idisk));

     // Bulge/bar fraction of 14 < I < 21
     // Assume same fMag for NSD as Bulge
     double fMagB =  fLF_detect(extI, Isst, Isen, 8);

     // Total
     double fMag = fMagD * f_disk + fMagB * (1 - f_disk);

     // Calc contribution of iDs to Nsall and tauall
     *Nsall  += (iDs>=2*narry21) ?            NDs[iDs]*fMag            
                                 :            NDs[iDs]*fMag*ks21[iDs/2];
     *tauall += (iDs>=2*narry21) ? tauDs[iDs]*NDs[iDs]*fMag 
                                 : tauDs[iDs]*NDs[iDs]*fMag*ks21[iDs/2];
     // printf ("%5d %.4f Ds[%02d]= %5.0f ls= %7.3f tauDs= %.4e sigN= %8.2f wtDB= %7.4f (m_idisk= %.1f f_disk= %.3f) f21D= %.4f - %.4f = %.4f f21B= %.4f - %.4f = %.4f sig21= %.4f tau*sig21= %.4e\n",idataAI,AIcoeff * (1 - exp(-Ds/hscale)),iDs,Ds,(Ds-Dsmin21)/dDs0,tauDs[iDs]*PI4GC2,NDs[iDs]*dDs0*STR2MIN2,wtDBs[iDs],m_idisk,f_disk,f21D,f14D,(f21D - f14D),f21B,f14B,(f21B - f14B),NDs[iDs]*f14_21*dDs0*STR2MIN2,tauDs[iDs]*NDs[iDs]*f14_21*PI4GC2*dDs0*STR2MIN2);
  }
  *Nsall *= dDs0;
  *tauall *= dDs0;
  // printf ("%5d %.4e\n",idataAI,tauall*PI4GC2);
  // Calc Nsall from Dlmin - Dsmin21
  double NDl  = 0;
  for (int iDl = 0; iDl < ncalcDl0; iDl++){    // loop for the next Dl integration
     int iDltmp = iDl - 2*narryDl + njiDl;        // iDltmp = njiDl - (nbunDl - njiDl)
     double Dl = (iDl>=2*narryDl) ? Dlmin + dDl0 * iDltmp 
               : (iDl   % 2 == 0) ? Dlmin + dDl0 * lsDl[iDl/2] : Dsmin21 - dDl0 * lsDl[iDl/2];

     // Extinction + Distance modulus at Dl
     double extI =  AI0 * (1 - exp(-Dl/hscale)) + 5 * log10(0.1*(Dl + 0.1));

     // DISK fraction of 14 < I < 21
     double m_idisk = (int) wtDBs[ncalc21+iDl];
     double f_disk  = 10*(wtDBs[ncalc21+iDl] - m_idisk);
     m_idisk *= 0.1;
     double wtD = m_idisk - (int) m_idisk; // e.g.) 4.2 - 4.0 = 0.2 
     double fMagD =  (1 - wtD) * fLF_detect(extI, Isst, Isen, (int)floor(m_idisk))
                   +      wtD  * fLF_detect(extI, Isst, Isen, (int) ceil(m_idisk));

     // Bulge/bar fraction of 14 < I < 21
     // Assume same fMag for NSD as Bulge
     double fMagB =  fLF_detect(extI, Isst, Isen, 8);

     // Total
     double fMag = fMagD * f_disk + fMagB * (1 - f_disk);

     // Calc contribution of i21 to Nsall and tauall
     NDl   += (iDl>=2*narryDl) ?            NDs[ncalc21+iDl]*fMag            
                               :            NDs[ncalc21+iDl]*fMag*ksDl[iDl/2];
     // printf ("%5d %.4f Dl[%02d]= %5.0f ls= %7.3f sigN= %8.2f wtDB= %7.4f (m_idisk= %.1f f_disk= %.3f) f21D= %.4f - %.4f = %.4f f21B= %.4f - %.4f = %.4f sigDl= %.4f\n",idataAI,AIcoeff * (1 - exp(-Dl/hscale)),iDl,Dl,(Dl-Dlmin)/dDl0,NDs[iDl+ncalc21]*dDl0*STR2MIN2,wtDBs[iDl+ncalc21],m_idisk,f_disk,f21D, f14D, (f21D - f14D),f21B, f14B,(f21B - f14B),NDs[iDl+ncalc21]*f14_21*dDl0*STR2MIN2);
  }
  NDl *= dDl0;
  // printf ("sig21= %.4f + %.4f = %.4f\n",Nsall*STR2MIN2,NDl*STR2MIN2,(Nsall+NDl)*STR2MIN2);
  *Nsall += NDl;
  *tauall = *tauall / *Nsall * PI4GC2;
  *Nsall *= STR2MIN2; // number / arcmin^2

} // End of optical depth calculation

//----------------
void store_NSDmoments(char *infile) // Read input_files/NSD_moments.dat
{
  // read moments of Sormani+21's NSD DF model 
  FILE *fp;
  char line[1000];
  char *words[100];
  if((fp=fopen(infile,"r"))==NULL){
     printf("can't open %s\n",infile);
     exit(1);
  }
  int iRz = 0;
  while (fgets(line,1000,fp) !=NULL){
     split((char*)" ", line, words);
     if (*words[0] == '#') continue;
     int iR = iRz % nRND;
     int iz = iRz / nRND;
     if (RstND + iR*dRND == 1000*atof(words[0]) && zstND + iz*dzND == 1000*atof(words[1])){
       logrhoNDs[iz][iR] = log10(atof(words[2])); // log [M_sun/pc^3]
       vphiNDs[iz][iR] = atof(words[3]); // vphi
       logsigvNDs[iz][iR][0] = log10(atof(words[4])); // sigphi
       logsigvNDs[iz][iR][1] = log10(atof(words[5])); // sigR
       logsigvNDs[iz][iR][2] = log10(atof(words[6])); // sigz
       corRzNDs[iz][iR] = atof(words[7]); // correlation coefficient between vR and vz
       // printf("iz=%d iR=%d %f %f %6.3f %5.1f\n", iz,iR,atof(words[1]),atof(words[0]),logrhoNDs[iz][iR], vphiNDs[iz][iR]);
     }else{
       printf("something goes wrong\n");
     }
     iRz++;
  } 
  fclose(fp);
}
//----------------
double like_obs(double mod, double obs, double err, double fe, int det, int UNIFORM){
  // return 0 when rejected, return 1 when accepted/not used, return 0 < f < 1 when fe > 0
  double Gamma_obs = 1;
  int sw = (det == 0) ? 1  // detection
         : (det == 1 && mod > obs) ? 1  // upper limit 
         : (det == 2 && mod < obs) ? 1  // lower limit 
         : 0;
  if(err >0 && sw == 1){
    if (UNIFORM == 1){
      Gamma_obs = (mod > obs-err && mod < obs+err) ? 1 : 0;
    }else{
      double chi, Plike, Plike0;
      chi = (mod - obs)/err;
      if (fe > 0){
        Plike0 = exp(-chi*chi*0.5);
        chi /= fe;
      }
      Plike = exp(-chi*chi*0.5);
      Gamma_obs = (Plike > ran1()) ? 1 : 0;
      if (fe > 0) Gamma_obs *= Plike0/Plike; // still 0 when not accepted
    }
  }
  return Gamma_obs;
}

//----------------
void store_IMF_nBs(int B, double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles, double M0, double M1, double M2, double M3, double Ml, double Mu, double alpha1, double alpha2, double alpha3, double alpha4, double alpha0){
  /* Store IMF with a broken-power law form.
   * Update normalize factors for the density distribution if B == 1 */
  double *PlogM_cum, *Mass, *PMlogM_cum, *PMlogM_cum_norm;
  Mass            = (double *)calloc(nm+1, sizeof(double *));
  PlogM_cum       = (double *)calloc(nm+1, sizeof(double *));
  PMlogM_cum      = (double *)calloc(nm+1, sizeof(double *));
  PMlogM_cum_norm = (double *)calloc(nm+1, sizeof(double *));
  logMst = log10(Ml);
	dlogM = (double) (log10(Mu)-logMst)/nm;
  for (int i=0; i<=nm; i++){
    double Mp  = i*dlogM + logMst;
    logMass[i] = Mp;
    Mass[i]  = pow(10, Mp);
    double alpha = (Mass[i] < M3) ? alpha4 : (Mass[i] < M2) ? alpha3 : (Mass[i] < M1) ? alpha2 : (Mass[i] < M0) ? alpha1 : alpha0;
    double temp00    = pow(M0,alpha0+1.);
    double temp01    = pow(M0,alpha1+1.);
    double temp11    = pow(M1,alpha1+1.);
    double temp12    = pow(M1,alpha2+1.);
    double temp22    = pow(M2,alpha2+1.);
    double temp23    = pow(M2,alpha3+1.);
    double temp33    = pow(M3,alpha3+1.);
    double temp34    = pow(M3,alpha4+1.);
    double dPlogM = 1;
    if (Mass[i] <M0) dPlogM=temp01 / temp00; //dM=MdlogM
    if (Mass[i] <M1) dPlogM=temp12 / temp11*dPlogM; //dM=MdlogM
    if (Mass[i] <M2) dPlogM=temp23 / temp22*dPlogM;
    if (Mass[i] <M3) dPlogM=temp34 / temp33*dPlogM;
    double templogMF = pow(Mass[i], alpha+1.);
    PlogM[i] = templogMF / dPlogM;
    if (i>=1) {
      PlogM_cum[i]  = 0.5*(PlogM[i]+PlogM[i-1])*dlogM                   + PlogM_cum[i-1];  // Mass function
      PMlogM_cum[i] = 0.5*(Mass[i]*PlogM[i]+Mass[i-1]*PlogM[i-1])*dlogM + PMlogM_cum[i-1]; // Mass spectrum
    } else {
      PlogM_cum[i] = 0.0;
      PMlogM_cum[i] = 0.0;
    }
  }
  // PlogM: Percentage of logM stars in total number of stars born, PMlogM: in total mass of stars born
  for(int i=0;i<=nm;i++){
    PlogM_cum_norm[i]= PlogM_cum[i]/PlogM_cum[nm];
    PMlogM_cum_norm[i]= PMlogM_cum[i]/PMlogM_cum[nm];
    PlogM[i] /= PlogM_cum[nm]; // for getcumu2xist
    int intp = PlogM_cum_norm[i]*20;
    if (imptiles[intp] == 0)  imptiles[intp] = (intp==0) ? 1 : i+0.5;
  }
  if (B == 0) return;

  // Calc average mass-loss for WDs
  double *ageMloss;
  ageMloss       = (double *)calloc(nm+1, sizeof(double *));
  double cumMwt = 0, cumWDwt = 0;
  void Mini2Mrem (double *pout, double M, int mean); 
  for (int i=nm;i>=0;i--){
    double pout[2] = {};
    double M = pow(10, logMass[i]);
    double wt = PlogM[i];
    Mini2Mrem(pout, M, 1);  // 0 : random
    double MWD = pout[0];
    cumMwt  += M * wt;
    cumWDwt += MWD * wt;
    ageMloss[i] = cumWDwt/cumMwt; 
  }
  // Read minimum died initial mass as a function of age
  char line[1000];
  char *words[100];
  FILE *fp;
  char file1[] = "input_files/Minidie.dat";
  double MRGstD[250], MRGenD[250], MRGstB[50], MRGenB[50], MRGstND[10], MRGenND[10];
  if((fp=fopen(file1,"r"))==NULL){
     printf("can't open %s\n",file1);
     exit(1);
  }
  nageD = 0, nageB = 0, nageND = 0;
  while (fgets(line,1000,fp) !=NULL){
     split((char*)" ", line, words);
     if (*words[0] == '#') continue;
     if (*words[0] == 'N'){
       agesND[nageND]    = atof(words[1]);
       MinidieND[nageND] = atof(words[2]);
       MRGstND[nageND] = atof(words[3]);
       MRGenND[nageND] = atof(words[4]);
       nageND++;
     }else if (*words[0] == 'B'){
       agesB[nageB]    = atof(words[1]);
       MinidieB[nageB] = atof(words[2]);
       MRGstB[nageB] = atof(words[3]);
       MRGenB[nageB] = atof(words[4]);
       nageB++;
     }else{
       agesD[nageD]    = atof(words[0]);
       MinidieD[nageD] = atof(words[1]);
       MRGstD[nageD] = atof(words[2]);
       MRGenD[nageD] = atof(words[3]);
       nageD++;
     }
  }
  fclose(fp);
  
  // for disks 
  double gamma = 1/tSFR;  // SFR timescale, 7 Gyr
  int agest = 1, ageen = 1000;
  int iages[7] = {15,100,200,300,500,700,1000};
  double wt_D[7] = {}, wtWD_D[7] = {}, sumM_D[7] = {}, sumMWD_D[7] = {}, sumstars_D[7] = {}, sumWDs_D[7] = {}, sumRGs_D[7] = {};
  for (int i=agest; i<=ageen; i++){
    int itmp = (i - agesD[0])/(agesD[1] - agesD[0]) + 0.5;
    if (itmp < 0) itmp = 0;
    double logMdie = log10(MinidieD[itmp]);
    double logMRG1 = log10(MRGstD[itmp]);
    double logMRG2 = log10(MRGenD[itmp]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG = PRG2 - PRG1; 
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1 - PM) * aveMloss;
    double PWD  = (1 - P);
    double wtSFR = exp(-gamma*(ageen-i)*0.01); // weight of this age
    P   *= wtSFR;
    PWD *= wtSFR;
    PM  *= wtSFR;
    PMWD *= wtSFR;
    PRG *= wtSFR;
    int idisk  = (i <= iages[0]) ? 0 
               : (i <= iages[1]) ? 1 
               : (i <= iages[2]) ? 2 
               : (i <= iages[3]) ? 3 
               : (i <= iages[4]) ? 4 
               : (i <= iages[5]) ? 5 
               : (i <= iages[6]) ? 6 : 0;
    wt_D[idisk]   += PM;
    wtWD_D[idisk] += PMWD;
    sumM_D[idisk] += PM*PMlogM_cum[nm];
    sumMWD_D[idisk]   += PMWD*PMlogM_cum[nm];
    sumstars_D[idisk] += P   *PlogM_cum[nm];
    sumWDs_D[idisk]   += PWD *PlogM_cum[nm];
    sumRGs_D[idisk]   += PRG *PlogM_cum[nm];
  }
  // Normalize
  double rho0thinMS = 0, rho0thinWD = 0, Sig2rho[8] = {}, aveMMS_D[8] = {}, aveMWD_D[8] = {}, nfracRG_D[8] = {}, aveM_D[8] = {};
  for (int i=0;i<8;i++){
    Sig2rho[i] = 0.5/zd[i]; // rho0/Sigma
    if (i < 7){
      int rd = (i==0) ? Rd[0] : Rd[1]; // because integrated mass depends on rd, SFR should be weight for the integrated mass  
      aveMMS_D[i] = sumM_D[i]/sumstars_D[i]; //  Msun/star for MainSequence
      aveMWD_D[i] = sumMWD_D[i]/sumWDs_D[i]; //  Msun/star for WhiteDwarf
      nfracRG_D[i]= sumRGs_D[i]/sumstars_D[i]; // RG to MS+RG ratio in number of stars
      aveM_D[i] = (sumM_D[i]+sumMWD_D[i])/(sumstars_D[i]+sumWDs_D[i]); // Msun/star for MS+WD
      // exp(-Rsun/rd)*wt[i]/rd is weight of rho at Sun position relative to the total mass wt[i] (but when ignoring hole)
      rho0thinMS += exp(-R0/rd)*wt_D[i]/rd * Sig2rho[i];
      rho0thinWD += exp(-R0/rd)*wtWD_D[i]/rd * Sig2rho[i];
    }
  }
  // double rhot0 = 0.042; // Msun/pc^3 @ z=0, rhot0 + rhoT0 = 0.042, (Bovy17: 0.042 +- 0.002 incl.BD)
  double rhoT0 = rhot0 * 0.04; //  Msun/pc^3, 4% of thin disk (Bland-Hawthorn & Gerhard (2016), f_rho = 4% +- 2%)
  for (int i=0;i<8;i++){
    int rd = (i == 0) ? Rd[0] : (i < 7) ? Rd[1] : (i == 7) ? Rd[2] : 0;
    if (i < 7){
      double norm = rhot0/rho0thinMS;
      double rhoMS  = norm * exp(-R0/rd) * wt_D[i]/rd * Sig2rho[i];
      double rhoWD  = norm * exp(-R0/rd) * wtWD_D[i]/rd * Sig2rho[i];
      rho0d[i] = rhoMS + rhoWD;
      n0MSd[i] = rhoMS/aveMMS_D[i];
      double n0WD = rhoWD/aveMWD_D[i];
      n0d[i]   = n0MSd[i] + n0WD;
      n0RGd[i] = n0MSd[i]*nfracRG_D[i];
    //   printf ("%d rho0= %.2e + %.2e = %.2e, n0= %.2e + %.2e = %.2e, n0RG= %.2e\n",i,rhoMS,rhoWD,rho0d[i],n0MSd[i],n0WD,n0d[i],n0RGd[i]);
    }else{  // Thick disk
      double logMdie = log10(MinidieD[nageD - 2]);
      double logMRG1 = log10(MRGstD[nageD - 2]);
      double logMRG2 = log10(MRGenD[nageD - 2]);
      double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
      double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
      double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
      double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
      double PRG = PRG2 - PRG1; 
      double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
      double PMWD = (1 - PM) * aveMloss;
      double PWD  = (1 - P);
      double aveMMS = PM   * PMlogM_cum[nm] / P   / PlogM_cum[nm]; // MSun/star for main sequence
      double aveMWD = PMWD * PMlogM_cum[nm] / PWD / PlogM_cum[nm]; // MSun/star for WD
      double aveM   = (PM*PMlogM_cum[nm]+PMWD*PMlogM_cum[nm])/(P*PlogM_cum[nm]+PWD*PlogM_cum[nm]);
      double norm = rhoT0/PM;
      double rhoMS = rhoT0;
      double rhoWD = norm * PMWD;
      rho0d[i] = rhoMS + rhoWD;
      n0MSd[i] = rhoMS/aveMMS;
      double n0WD = rhoWD/aveMWD;
      n0d[i]   = n0MSd[i] + n0WD;
      n0RGd[i] = n0MSd[i]* PRG/P;
      // printf ("%d rho0= %.2e + %.2e = %.2e, n0= %.2e + %.2e = %.2e, n0RG= %.2e\n",i,rhoMS,rhoWD,rho0d[i],n0MSd[i],n0WD,n0d[i],n0RGd[i]);
    }
  }
  // for Bar 
  // Use 9+-1 Gyr to calculate the conversion factors (e.g., for total mass -> MS).
  // This is because the K21 fit was done with this assamption.
  // In the calculation of magnitude or PDMF, only the isochrone with 9Gyr is used, though.
  // So, a small discrepancy exists in the total mass to MS mass ratio between fb_MS value and MC simulation
  double wt_B = 0, wtWD_B = 0, sumM_B = 0, sumMWD_B = 0, sumstars_B = 0, sumWDs_B = 0, sumRGs_B = 0;
  for (int i= 0; i< nageB; i++){
    double tau = 0.01*agesB[i];
    double wtSFR = (tau - mageB)/sageB;
    wtSFR = exp(-0.5*wtSFR*wtSFR);
    double logMdie = log10(MinidieB[i]);
    double logMRG1 = log10(MRGstB[i]);
    double logMRG2 = log10(MRGenB[i]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG = PRG2 - PRG1; 
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1 - PM) * aveMloss;
    double PWD  = (1 - P);
    P   *= wtSFR;
    PWD *= wtSFR;
    PM  *= wtSFR;
    PMWD *= wtSFR;
    PRG *= wtSFR;
    wt_B   += PM;
    wtWD_B += PMWD;
    sumM_B += PM*PMlogM_cum[nm];
    sumMWD_B   += PMWD*PMlogM_cum[nm];
    sumstars_B += P   *PlogM_cum[nm];
    sumWDs_B   += PWD *PlogM_cum[nm];
    sumRGs_B   += PRG *PlogM_cum[nm];
  }
  double aveMMS = sumM_B/sumstars_B;
  double aveMWD = sumMWD_B/sumWDs_B;
  double aveM   = (sumM_B+sumMWD_B)/(sumstars_B+sumWDs_B);
  m2nb_MS  = 1/aveMMS;
  m2nb_WD  = 1/aveMWD;
  nMS2nRGb = sumRGs_B/sumstars_B; // RG to MS+RG ratio in number of stars
  fb_MS    = wt_B/(wt_B+wtWD_B);

  // for NSD
  double wt_ND = 0, wtWD_ND = 0, sumM_ND = 0, sumMWD_ND = 0, sumstars_ND = 0, sumWDs_ND = 0, sumRGs_ND = 0;
  // As of 20220207, nageND = 1.
  for (int i= 0; i< nageND; i++){
    double tau = 0.01*agesND[i];
    double wtSFR = (tau - mageND)/sageND;
    wtSFR = exp(-0.5*wtSFR*wtSFR);
    double logMdie = log10(MinidieND[i]);
    double logMRG1 = log10(MRGstND[i]);
    double logMRG2 = log10(MRGenND[i]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG = PRG2 - PRG1; 
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1 - PM) * aveMloss;
    double PWD  = (1 - P);
    P   *= wtSFR;
    PWD *= wtSFR;
    PM  *= wtSFR;
    PMWD *= wtSFR;
    PRG *= wtSFR;
    wt_ND   += PM;
    wtWD_ND += PMWD;
    sumM_ND += PM*PMlogM_cum[nm];
    sumMWD_ND   += PMWD*PMlogM_cum[nm];
    sumstars_ND += P   *PlogM_cum[nm];
    sumWDs_ND   += PWD *PlogM_cum[nm];
    sumRGs_ND   += PRG *PlogM_cum[nm];
  }
  aveMMS = sumM_ND/sumstars_ND;
  aveMWD = sumMWD_ND/sumWDs_ND;
  aveM   = (sumM_ND+sumMWD_ND)/(sumstars_ND+sumWDs_ND);
  m2nND_MS  = 1/aveMMS;
  m2nND_WD  = 1/aveMWD;
  nMS2nRGND = sumRGs_ND/sumstars_ND; // RG to MS+RG ratio in number of stars
  fND_MS    = wt_ND/(wt_ND+wtWD_ND);
  free(Mass          );
  free(PlogM_cum     );
  free(PMlogM_cum    );
  free(PMlogM_cum_norm);
}

//----------------
void Mini2Mrem (double *pout, double Mini, int mean) {  // mean = 1: give mean, 0: give random
  /* Return remnant mass for a given initial mass following the initial-final mass relation by Lam et al. 2020, ApJ, 889, 31 */
  double Mrem, fREM; 
  // Below is from Table 1 of Lam et al. 2020, ApJ, 889, 31
  double PNS = (Mini < MiniWDmax) ? 0  // 100% WD
             : (Mini < 15.0) ? 1  // 100% NS
             : (Mini < 17.8) ? 0.679
             : (Mini < 18.5) ? 0.833
             : (Mini < 21.7) ? 0.500
             : (Mini < 25.2) ? 0  // 100% BH
             : (Mini < 27.5) ? 0.652 
             : (Mini < 60.0) ? 0  // 100% BH
             : 0.4;
  // IFMRs for NS and BH are from Appendix C of Lam et al. 2020, ApJ, 889, 31 or from Raithel+18
  if (Mini < MiniWDmax){
     Mrem = 0.109 * Mini + 0.394; // IFMR from Kalirai+08
     fREM = 1; // WD
  }else{ 
     // NS (Eqs.(11)-(16) of Raithel+18)
     double MNS;
     do {
       MNS = (Mini < 13.0) ? 2.24 + 0.508 *(Mini - 14.75) 
                                  + 0.125 *(Mini - 14.75)*(Mini - 14.75) 
                                  + 0.011 *(Mini - 14.75)*(Mini - 14.75)*(Mini - 14.75)
           : (Mini < 15.0) ?  0.123 + 0.112 * Mini
           : (Mini < 17.8) ?  0.996 + 0.0384* Mini
           : (Mini < 18.5) ? -0.020 + 0.10  * Mini
           : (Mini < 21.7 && mean == 0) ? 1.60 + 0.158*gasdev()
           : (Mini < 21.7 && mean == 1) ? 1.60 
           : (Mini < 27.5) ?  3232.29 - 409.429*(Mini - 2.619) 
                                      + 17.2867*(Mini - 2.619)*(Mini - 2.619) 
                                      - 0.24315*(Mini - 2.619)*(Mini - 2.619)*(Mini - 2.619)
           : (mean == 0) ? 1.78 + 0.02*gasdev()
           : 1.78;
           // print "Mini=Mini, MNS = MNS\n";
     }while(PNS > 0 && (MNS < MNSMIN || MNS > MNSMAX));
     // BH
     double Mcore = (Mini < 42.21) ? -2.049 + 0.4140 * Mini
                 : 5.697 + 7.8598 * 1e+8 * pow(Mini, -4.858);
     double Mall = 15.52 - 0.3294*(Mini - 25.97) 
                       - 0.02121*(Mini - 25.97)*(Mini - 25.97) 
                      + 0.003120*(Mini - 25.97)*(Mini - 25.97)*(Mini - 25.97);
     double fej = (Mini < 42.21) ? 0.9 : 1.0;
     double MBH = fej*Mcore + (1-fej)*Mall;
     // print "Mini=Mini, MBH = MBH\n";

     // Mean or Rand
     if (mean == 1){
       Mrem = PNS*MNS + (1-PNS)*MBH;
       fREM = PNS*2 + (1-PNS)*3;
     }else{
       double ran = ran1();
       Mrem = (ran < PNS) ? MNS : MBH;
       fREM = (ran < PNS) ?    2 :    3;
     }
  }
  pout[0] = Mrem;
  pout[1] = fREM;
}
//----------------
double fLF_detect(double extI, double Imin, double Imax, int idisk){
  double imaxd = (Imax - extI - MIs[0])/dILF;
  double imind = (Imin - extI - MIs[0])/dILF;
  if (imaxd < 0)      imaxd = 0;
  if (imaxd > nMIs-1) imaxd = nMIs - 1;
  if (imind < 0)      imind = 0;
  if (imind > nMIs-1) imind = nMIs - 1;
  int imax = imaxd;
  int imin = imind;
  double fmax = CumuN_MIs[idisk][imax+1]*(imaxd-imax) 
              + CumuN_MIs[idisk][imax]  *(1 - (imaxd-imax)); 
  double fmin = CumuN_MIs[idisk][imin+1]*(imind-imin) 
              + CumuN_MIs[idisk][imin]  *(1 - (imind-imin)); 
  return (fmax - fmin);
}
//----------------
double fIVI_detect(double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk){
  double imaxd = (Imax - extI - MIs[0])/dILF;
  double imind = (Imin - extI - MIs[0])/dILF;
  if (imaxd < 0)      imaxd = 0;
  if (imaxd > nMIs-1) imaxd = nMIs - 1;
  if (imind < 0)      imind = 0;
  if (imind > nMIs-1) imind = nMIs - 1;
  int imax = imaxd;
  int imin = imind;
  double jmaxd = (VImax - extVI - VIs[0])/dVILF;
  double jmind = (VImin - extVI - VIs[0])/dVILF;
  if (jmaxd < 0)      jmaxd = 0;
  if (jmaxd > nVIs-1) jmaxd = nVIs - 1;
  if (jmind < 0)      jmind = 0;
  if (jmind > nVIs-1) jmind = nVIs - 1;
  int jmax = jmaxd;
  int jmin = jmind;
  double fIVI = 0;
  for (int j= jmin; j <= jmax; j++){
  for (int i= imin; i <= imax; i++){
    fIVI += f_VI_Is[idisk][j][i];
    // printf("VI= %.3f MI= %.3f f= %.5e sumf= %.5e\n",VIs[j],MIs[i],f_VI_Is[idisk][j][i],fIVI);
  }}
  return fIVI;
}
//----------------
void store_cumuP_Shu(char *infile) // calculate cumu prob dist of fg = Rg/R following Shu DF
{
  // read circular velocity
  FILE *fp;
  char line[1000];
  char *words[100];
  if (nVcs == 0){
    if((fp=fopen(infile,"r"))==NULL){
       printf("can't open %s\n",infile);
       exit(1);
    }
    nVcs = 0;
    while (fgets(line,1000,fp) !=NULL){
       split((char*)" ", line, words);
       if (*words[0] == '#') continue;
       Rcs[nVcs]  = 1000*atof(words[0]); // kpc -> pc
       Vcs[nVcs] =      atof(words[1]); // km/sec
       nVcs++;
    } 
    fclose(fp);
  }
  // Store CPD of fg following Shu DF
  // v[iz][iR][idisk]
  double getx2y(int n, double *x, double *y, double xin);
  double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd);
  void get_PRRGmax2(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
  for (int z = zstShu; z <= zenShu; z+=dzShu){
    int iz = (z - zstShu)/dzShu;
    double facVcz = 1 + 0.0374*pow(0.001*fabs(z), 1.34); // Eq. (22) of Sharma et al. 2014, ApJ, 793, 51
    for (int R = RstShu; R <= RenShu; R+=dRShu){
      int iR = (R - RstShu)/dRShu;
      double vcR  = getx2y(nVcs, Rcs, Vcs, R);
      for (int idisk=0; idisk<8; idisk++){
        double tau = medtauds[idisk];
        double hsigU = (idisk < 7) ? hsigUt : hsigUT;
        int    rd = (idisk == 0) ? Rd[0] : (idisk <  7) ? Rd[1] : Rd[2];
        double sigU0 = (idisk < 7) ? sigU10d * pow((tau+0.01)/10.01, betaU) : sigU0td;
        double Rgmin = R0 - hsigU*log(vcR/sigU0); // which gives c = 0.5 if vcR = vcRg
        if (Rgmin > R) Rgmin = R0 - hsigU*log(240.0/sigU0); // vcmax = 240
        double fgmin0 = Rgmin/R;
        double fg1 = (fgmin0 > 1.5) ? fgmin0 : 1; // initial value of Newton method in get_PRRGmax
        double pout[4] = {};
        get_PRRGmax2(pout, R, z, fg1, sigU0, hsigU, rd);
        double   Pmax = pout[0];
        double  fgmin = pout[1];
        double  fgmax = pout[2];
        double    fgc = pout[3];
        if ((fgmin > 1 && R > 1000) || Pmax == 0) 
          printf ("# PERROR!! get_PRRGmax2(pout, %5d, %4d, %.3f, %.2f, %.2f, %d)\n",R, z, fg1, sigU0, hsigU, rd);
        // if (fgmin < 0.1 && fgc > 0.5) fgmin = 0.1;
        int swerror = ((fgmin > 1 && R > 1000) || Pmax == 0) ? 1 : 0;
        double fg   = fgmin;
        double dfg0 = (fgc - fgmin)*0.025; // divided by 40
        int ifg = 0; 
        double dfg =0;
        while(fg <= fgmax){
          fgsShu[iz][iR][idisk][ifg] = fg;
          double PRRg = calc_PRRg(R,z,fg,sigU0,hsigU,rd);
          PRRgShus[iz][iR][idisk][ifg] = PRRg;
          cumu_PRRgs[iz][iR][idisk][ifg] = (ifg==0) ? 0 : cumu_PRRgs[iz][iR][idisk][ifg-1] + 0.5*(PRRgShus[iz][iR][idisk][ifg-1] + PRRgShus[iz][iR][idisk][ifg])*dfg;
          dfg = (PRRg/Pmax < 0.05) ? 4*dfg0 : (PRRg/Pmax < 0.25 || PRRg/Pmax > 0.7) ? dfg0 : 2*dfg0;
          //  idfg = (abs(fgc-fg) <= 0.10) ? 0.02 : 0.06;
          //  printf "%2d (%.3f)  %.4f %.5e %.5e\n",ifg,fgmin,fg,PRRg,cumu_PRRgs[iz][iR][idisk][ifg]; 
          ifg++;
          fg = fg + dfg;
        }
        n_fgsShu[iz][iR][idisk] = ifg;
        // normalize and store percentiles
        double norm = cumu_PRRgs[iz][iR][idisk][ifg-1];
        for (int ktmp=0; ktmp<ifg;ktmp++){
          PRRgShus[iz][iR][idisk][ktmp]   /= norm;
          cumu_PRRgs[iz][iR][idisk][ktmp] /= norm;
          int intp = cumu_PRRgs[iz][iR][idisk][ktmp]*20;
          if (kptiles[iz][iR][idisk][intp]==0) kptiles[iz][iR][idisk][intp] = (intp==0) ? 1 : ktmp+0.5;
          // printf("(%4d-%4d-%d) ktmp= %3d (< %3d), fg= %.3f PRRg= %.4e (f= %.4f) cumu_PRRg= %.4e intp= %2d, kptile[intp]= %2d\n",z,R,idisk,ktmp,ifg,fgsShu[iz][iR][idisk][ktmp],PRRgShus[iz][iR][idisk][ktmp],PRRgShus[iz][iR][idisk][ktmp]/(Pmax/norm),cumu_PRRgs[iz][iR][idisk][ktmp], intp, kptiles[iz][iR][idisk][intp]);
        }
        if (swerror == 1) 
           printf("# i=%d, tau=%5.2f fg= %7.4f - %7.4f, fgc= %6.4f Pmax= %.3e\n",idisk,tau,fgmin,fgmax,fgc,Pmax);
      }
    }
  }
}
//---- calc Pmax, fgmin, fgmax, fgc -------
void get_PRRGmax2(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd){
  if (fg1 < 1) fg1 = 1;
  double dfg = 0.001;
  double fgc = 1e+3, Pmax = 1e-200, dPdfgc = 0;
  double Ptmp = 0;
  double fg, fg2, fg3, fg4, dPdfg1, dPdfg2, d2Pdfg, dPdfg3, dPdfg4, d2Pdfg2, P1, P2, P3, P4;
  double jj;
  int    nj = 0, ntry = 0, sw = 0;
  double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd);
  void calc_dpdfg(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
  if (hsigU/rd/sigU0 < 0.1){  // eg hsigU/rd = 4 && sigU0 > 40 or hsigU/rd = 3.5 && sigU0 > 35
    for (fg=0.15;fg<1.0;fg+=0.05){ // determine fgmin & fgmax from Newton's method (search P = 0)
      P1 = calc_PRRg(R,z,fg,sigU0,hsigU,rd);
      if (P1 > Ptmp){
        Ptmp = P1;
        fg1 = fg;
      }
    }
  }
  while(1){
    int ncalc = 0;
    double pout1[2] = {}, pout2[2]= {};
    for (int j=0;j<3;j++){ // Find Pmax by Newton's method (search dP/dfg = 0)
      fg2 = fg1 + dfg;
      calc_dpdfg(pout1,R,z,fg1,sigU0,hsigU,rd); // d(PRRg)/d(fg)
      calc_dpdfg(pout2,R,z,fg2,sigU0,hsigU,rd);
      dPdfg1 = pout1[0], dPdfg2 = pout2[0];
      P1     = pout1[1], P2     = pout2[1];
      d2Pdfg = (dPdfg2-dPdfg1)/dfg; // d2(PRRg)/d(fg)
      if (P1 > Pmax){ fgc    = fg1    ;
                      dPdfgc = dPdfg1 ;
                      Pmax   = P1     ;}
      // printf "# R=%5d, now(j, ncalc, fg1, dPdfg1, dPdfg2, d2Pdfg, P)= (j, %3d, %.3f, %8.1e, %8.1e, %8.1e, %.4e), best(fgc, dPdfgc, Pmax)= (%.3f, %8.1e, %.4e)\n",R,ncalc,fg1,dPdfg1,dPdfg2,d2Pdfg,P1,fgc,dPdfgc,Pmax;
      ncalc++;
      if (ncalc > 15){
        if (nj > 0){
          break;
        }else if(ntry < 2){
          if (fgc > 900) fgc = (ntry == 0) ? fg1 : 0.9;
          fg1 = (ntry == 0) ? fgc - 0.4 : fgc + 0.4;
          if (fg1 < 0) fg1 = 0.2*ran1();
          ncalc = 0; 
          j = -1;
          ntry++;
          continue;
        }else{
          // printf ("# break!!\n");
          break;
        }
      }
      if (j==2 && fabs(dPdfgc/Pmax) > 0.1){
        nj++;
        fg1 = (dPdfgc > 0) ? fgc + 0.05/nj*ran1() : fgc - 0.05/nj*ran1();
        j = -1;
        continue;
      }
      if (dPdfg1 == 0){ // too left or too right
        jj = (dPdfgc == 0) ? 0.5 : 0.2*ran1();
        fg1 = (fg1 < fgc) ? fg1 + jj : fg1 - jj;
        j=-1;
        continue;
      }
      if (d2Pdfg > 0 && dPdfg1 < 0){ // to confirm too right or marume
        fg3 = fg2 + 0.04; // 0.04 ha tekitou
        fg4 = fg3 + dfg;
        calc_dpdfg(pout1,R,z,fg3,sigU0,hsigU,rd);
        calc_dpdfg(pout2,R,z,fg4,sigU0,hsigU,rd);
        dPdfg3 = pout1[0], dPdfg4 = pout2[0];
        P3     = pout1[1], P4     = pout2[1];
        d2Pdfg2 = (dPdfg4-dPdfg3)/dfg;
        if (d2Pdfg2 > 0 || dPdfg3 == 0){ // too right
          fg1 -= (0.02 + 0.10*ran1());
          j=-1;
          continue;
        }else{
          d2Pdfg = d2Pdfg2; // d2Pdfg > 0 was marume and replaced with d2Pdfg2 
        }
      }
      if (d2Pdfg > 0 && dPdfg1 > 0){ // to confirm too left or marume
        fg3 = fg1 - 0.04; // 0.04 ha tekitou
        fg4 = fg3 + dfg;
        calc_dpdfg(pout1,R,z,fg3,sigU0,hsigU,rd);
        calc_dpdfg(pout2,R,z,fg4,sigU0,hsigU,rd);
        dPdfg3 = pout1[0], dPdfg4 = pout2[0];
        P3     = pout1[1], P4     = pout2[1];
        d2Pdfg2 = (dPdfg4-dPdfg3)/dfg;
        if (d2Pdfg2 > 0 || dPdfg3 == 0){ // too left
          fg1 += (0.02 + 0.10*ran1());
          j=-1;
          continue;
        }else{
          d2Pdfg = d2Pdfg2; // d2Pdfg > 0 was marume and replaced with d2Pdfg2 
        }
      }
      // printf "# R=%5d, fg1(j)= %.3f, dPdfg1= %8.1e, dPdfg2= %8.1e, d2Pdfg= %8.1e, P=%.4e",R,fg1,dPdfg1,dPdfg2,d2Pdfg,P1;
      if (d2Pdfg != 0) fg1 = fg1 - dPdfg1/d2Pdfg;
      if (fg1 < 0) fg1  = 0.1;
      if (fabs(dPdfg1/d2Pdfg) > 0.5){
        jj = (dPdfgc > 0) ? 0.10 : -0.10;
        fg1 = fgc + jj*ran1();
        j = -1;
        continue;
      }
      // printf " fgc=fgc, Pmax=Pmax, nextfg1= %.3f\n",fg1;
    }
    ncalc = 0; sw = 0;
    for (fg1=fgc-0.2;fg1>0.1;fg1-=0.2){ // determine fgmin & fgmax from Newton's method (search P = 0)
      P1 = calc_PRRg(R,z,fg1,sigU0,hsigU,rd);
      ncalc++;
      if (P1 > Pmax*1.05){
        // printf "# P1 (P1 @ fg1) > Pmax (Pmax @ fgc)!! Calc again!\n";
        Pmax = P1;
        fgc = fg1;
        sw = 1;
      }
      // printf "fg1 -> %.2e (%.2e)\n",P1,P1/Pmax;
      if (P1/Pmax < 1e-02) break;
    }
    if (sw == 1) fg1 = fgc;
    // print "# next w/ fg1 = fg1 = fgc\n" if sw ==1;
    if (sw == 1) continue;
    sw = 0;
    for (fg2=fgc+0.2;fg2<4.0;fg2+=0.2){ // determine fgmin & fgmax from Newton's method (search P = 0)
      P2 = calc_PRRg(R,z,fg2,sigU0,hsigU,rd);
      ncalc++;
      if (P2 > Pmax*1.05){
        Pmax = P2;
        fgc = fg2;
        sw = 1;
        break;
      }
      // printf "fg2 -> %.2e (%.2e)\n",P2,P2/Pmax;
      if (P2/Pmax < 1e-02) break;
    }
    if (sw == 1) fg1 = fgc;
    // print "# next w/ fg1 = fg1 = fgc\n" if sw ==1;
    if (sw == 1) continue;
    if (fg1 < 0) fg1 = 0.1;
    // printf "#fg1= fg1, fg2=fg2, ncalc=ncalc, fP1=%.5f, fP2= %.5f\n",P1/Pmax,P2/Pmax;
    pout[0] = Pmax;
    pout[1] =  fg1;
    pout[2] =  fg2;
    pout[3] =  fgc;
    break;
  }
}
//---- calc P(Rg|R) following Shu distribution ( Eq.(16) of BG16 ) -------
void calc_dpdfg(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd){
  double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd);
  double dfg = 0.001;
  double fg2 = fg1 + dfg;
  double PRRg1 = calc_PRRg(R,z,fg1,sigU0,hsigU,rd);
  double PRRg2 = calc_PRRg(R,z,fg2,sigU0,hsigU,rd);
  double dPdfg = (PRRg2-PRRg1)/dfg;
  if (PRRg2 <= 0 || PRRg1 <= 0){
    dPdfg = PRRg1 = 0;  
  }
  // print "fg1=fg1, vc1=vc1, a01=a01, a1= a1, P=PRRg1\n";
  // print "fg2=fg2, vc2=vc2, a02=a02, a2= a2, P=PRRg2\n";
  pout[0] = dPdfg;
  pout[1] = PRRg1;
}
//----------------
void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD) //
{
  void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);
  double getcumu2xist (int n, double *x, double *F, double *f, double Freq, int ist, int inv);
  double getx2y(int n, double *x, double *y, double xin);
  double xyz[3]={};
  Dlb2xyz(D, lD, bD, R0, xyz);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  double R = sqrt(x*x + y*y);
  double vx = 0, vy = 0, vz = 0;
  if (i < 8){
    double sigW0 = (i < 7) ? sigW10d * pow((tau+0.01)/10.01, betaW) : sigW0td;
    double sigU0 = (i < 7) ? sigU10d * pow((tau+0.01)/10.01, betaU) : sigU0td;
    double hsigW = (i < 7) ? hsigWt : hsigWT;
    double hsigU = (i < 7) ? hsigUt : hsigUT;
    double sigW  = sigW0*exp(-(R - R0)/hsigW);
    double sigU  = sigU0*exp(-(R - R0)/hsigU);
    int iz = (fabs(z) - zstShu)/dzShu;
    int iR = (R > RstShu) ? (R - RstShu)/dRShu : 0; // R = RstShu if R < RstShu
    do{
      double ran = ran1();
      int inttmp = ran*20;
      int kst1 = 1, kst2 = 1, kst3 = 1, kst4 = 1; // to avoid bug when inttmp = 0
      for (int itmp = inttmp; itmp > 0; itmp--){
        if (kst1 == 1) kst1 = kptiles[iz][iR][i][itmp];
        if (kst1 > 0 && kst2 > 0 && kst3 > 0 && kst4 > 0) break;
      }
      double fg1= getcumu2xist(n_fgsShu[iz][iR][i]    , fgsShu[iz][iR][i]    ,cumu_PRRgs[iz][iR][i]    ,PRRgShus[iz][iR][i]    ,ran,kst1,0);
      double fg = fg1;
      double Rg = fg*R;
      double vc = getx2y(nVcs, Rcs, Vcs, Rg) / (1 + 0.0374*pow(0.001*fabs(z), 1.34));
      double vphi = vc*fg;
      double vR =    0 + gasdev()*sigU; // radial velocity
      vx = -vphi * y/R + vR * x/R; // x/R = cosphi, y/R = sinphi
      vy =  vphi * x/R + vR * y/R; // x/R = cosphi, y/R = sinphi
      vz =    0 + gasdev()*sigW; // vertical velocity
    }while (vx*vx + vy*vy + vz*vz > vescd*vescd);
  }else if (i == 9 && ND == 3){ // NSD (when ND == 3)
    if (R > RenND || fabs(z) > zenND){
      printf("ERROR: NSD comp exists where it must not exist. (R,z)= (%f, %f)!!\n",R,z);
      exit(1);
    }
    // Bilinear interpolation of Sormani+21's NSD DF moments
    double as[4] = {}; // coeffs for interpolation
    double m_vphi = 0, logsigphi = 0, logsigR = 0, logsigz = 0, corRz = 0;
    interp_xy_coeff(nzND, nRND, as, zstND, RstND, dzND, dRND, fabs(z), R);
    int iz0  = (fabs(z) - zstND)/dzND;
    int iR0  = (R - RstND)/dRND;
    for (int j = 0; j < 4; j++){
      int iz = (j == 0 || j == 2) ? iz0 : iz0 + 1;
      int iR = (j == 0 || j == 1) ? iR0 : iR0 + 1;
      if (as[j] > 0){
        m_vphi    += as[j]*vphiNDs[iz][iR];
        logsigphi += as[j]*logsigvNDs[iz][iR][0];
        logsigR   += as[j]*logsigvNDs[iz][iR][1];
        logsigz   += as[j]*logsigvNDs[iz][iR][2];
        corRz     += as[j]*corRzNDs[iz][iR];
        // printf ("%d %d %d %f %f %f %f %f\n",j,iz,iR,m_vphi,logsigphi,logsigR,logsigz,corRz);
      }
    }
    // Random velocity with correlation coeff between vR and vz
    // ref: https://www.sas.com/offices/asiapacific/japan/service/technical/faq/list/body/stat034.html
    double sigphi = pow(10.0, logsigphi);
    double sigR   = pow(10.0, logsigR);
    double sigz   = pow(10.0, logsigz);
    double facR   = sigz/sigR * corRz;
    double sigz_R = sigz*sqrt(1 - corRz*corRz);
    do{
      double vphi = m_vphi + gasdev()*sigphi; // Assume vphi distribution is symmetrical (which is not true)
      double vR = gasdev()*sigR;
      vx = -vphi * y/R + vR * x/R; // x/R = cosphi, y/R = sinphi
      vy =  vphi * x/R + vR * y/R; // x/R = cosphi, y/R = sinphi
      vz =  facR * vR  + gasdev()*sigz_R; // random w/ correlation coeff
      // printf("%f %f %f %f %f %f %f %f %f %f\n",R,z,vphi,vR,vz,corRz,m_vphi,sigphi,sigR,sigz);
    }while (vx*vx + vy*vy + vz*vz > vescb*vescb);
  }else{ // bar & NSD (when ND <= 2)
    double vrot = 0.001 * Omega_p * R; // km/s/kpc -> km/s/pc
    double xb =  x * costheta + y * sintheta;
    double yb = -x * sintheta + y * costheta;
    double zb =  z;                          
    double sigvbs[3] = {}, sigx, sigy, sigz;
    void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
    calc_sigvb(xb, yb, zb, sigvbs);
    sigx = sqrt(sigvbs[0]*sigvbs[0] * costheta*costheta + sigvbs[1]*sigvbs[1] * sintheta*sintheta);
    sigy = sqrt(sigvbs[0]*sigvbs[0] * sintheta*sintheta + sigvbs[1]*sigvbs[1] * costheta*costheta);
    sigz = sigvbs[2];
    double avevxb   = (yb > 0) ? -vx_str : vx_str;
    if (y0_str > 0){
      double tmpyn = fabs(yb/y0_str);
      avevxb  *=  (1 - exp(-tmpyn*tmpyn));
    }
    do{
      vx = - vrot * y/R + avevxb * costheta + sigx * gasdev();
      vy =   vrot * x/R + avevxb * sintheta + sigy * gasdev();
      vz =                                    sigz * gasdev();
    }while (vx*vx + vy*vy + vz*vz > vescb*vescb);
  }
  vxyz[0] = vx;
  vxyz[1] = vy;
  vxyz[2] = vz;
}
//---------------
void getaproj(double *pout, double M1, double M2, int coeff)  { // pick up aproj  
   double Mprim = (M2>M1) ? M2 : M1;
   double meanloga  = 0.57 + 1.02 * Mprim; // Table 2 of Koshimoto+20, AJ, 159, 268
   if (meanloga > MAXMEANLOGA) meanloga = MAXMEANLOGA; // avoid too small meanmaloga 
   if (meanloga < MINMEANLOGA) meanloga = MINMEANLOGA; // avoid too small meanmaloga 
   double sigmaloga = 1.61 + 1.15 * log(Mprim)/log(10); // Table 2 of Koshimoto+20, AJ, 159, 268
   if (sigmaloga > MAXSIGLOGA) sigmaloga = MAXSIGLOGA; // avoid too small sigmaloga 
   if (sigmaloga < MINSIGLOGA) sigmaloga = MINSIGLOGA; // avoid too small sigmaloga 
   double ran = coeff * fabs(gasdev());
   double loga = meanloga + ran * sigmaloga;
   // logaproj = loga - 0.133; // 0.133 = <log(a/aproj)>
   double a = pow(10.0, loga);
   ran = ran1();
   double aproj = sqrt(1 - ran*ran) * a; //probability of rproj/a from Gould&Loeb 1992
   pout[0] = loga;
   pout[1] = aproj;
}

double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv){ 
  // for cumulative distribution (assuming linear interpolation for f(x) when cumu = F = int f(x))
  double Fmax = F[n-1];
  double Fmin = F[0];
  if (Fmin > Freq) return 0;
  if (Fmax < Freq) return 0;
  if (ist < 1) ist = 1;
  if (inv==0){
    for(int i=ist;i<n;i++){
       if (F[i] <= Freq && F[i-1] >Freq || F[i] >=Freq && F[i-1] < Freq){
          double a = 0.5*(f[i]-f[i-1])/(x[i]-x[i-1]);
          double b = f[i-1] - 2*a*x[i-1];
          double c = a*x[i-1]*x[i-1] - f[i-1]*x[i-1] + F[i-1] - Freq;
          double xreq = (a != 0) ? (-b + sqrt(b*b - 4*a*c)) * 0.5/a  // root of ax^2 +bx + c
                                 : (x[i]-x[i-1])/(F[i]-F[i-1])*(Freq-F[i-1]) + F[i-1];
          return xreq;
       }
    }
  }else{
    for(int i=ist;i>0;i--){
       if (F[i] <= Freq && F[i-1] >Freq || F[i] >=Freq && F[i-1] < Freq){
          double a = 0.5*(f[i]-f[i-1])/(x[i]-x[i-1]);
          double b = f[i-1] - 2*a*x[i-1];
          double c = a*x[i-1]*x[i-1] - f[i-1]*x[i-1] + F[i-1] - Freq;
          double xreq = (a != 0) ? (-b + sqrt(b*b - 4*a*c)) * 0.5/a  // root of ax^2 +bx + c
                                 : (x[i]-x[i-1])/(F[i]-F[i-1])*(Freq-F[i-1]) + F[i-1];
          return xreq;
       }
    }
  }
  return 0;
}
//---------------
int read_MLemp(char *infile, double *M_emps, double **Mag_emps) 
{
   FILE *fp;
   char line[1000];
   char *words[100];
   if((fp=fopen(infile,"r"))==NULL){
      printf("can't open %s\n",infile);
      exit(1);
   }
   int i=0;
   while (fgets(line,1000,fp) !=NULL){
      int nwords = split((char*)" ", line, words);
      if (*words[0] == '#') continue;
      M_emps[i] = atof(words[0]);
      // printf ("# %2d %.3f",i,M_emps[i]);
      for (int j=1; j<nwords; j++){ 
        Mag_emps[j-1][i] = atof(words[j]);
        // printf (" %6.3f",Mag_emps[j-1][i]);
      }
      // printf ("\n");
      i++;
   }
   fclose(fp);
   return i;
}
//---------------
int make_LFs(double *MIs, double **CumuN_MIs, double *logMass, double *PlogM_cum_norm) 
{
   char *infile;
   FILE *fp;
   char line[1000];
   char *words[100];
   int i=0;
   infile = (char*)"input_files/NbleNall_bin.dat";
   if((fp=fopen(infile,"r"))==NULL){
      printf("can't open %s\n",infile);
      exit(1);
   }
   while (fgets(line,1000,fp) !=NULL){
      int nwords = split((char*)" ", line, words);
      if (*words[0] == '#') continue;
      MIs[i] = atof(words[0]);
      i++;
   }
   fclose(fp);

   // Make LFs
   const char *file1;
   int nage, narry;
   for (int icomp=0; icomp<ncomp; icomp++){
     file1 = (icomp == 0) ? "input_files/iso-thin1.dat" :
             (icomp == 1) ? "input_files/iso-thin2.dat" :
             (icomp == 2) ? "input_files/iso-thin3.dat" :
             (icomp == 3) ? "input_files/iso-thin4.dat" :
             (icomp == 4) ? "input_files/iso-thin5.dat" :
             (icomp == 5) ? "input_files/iso-thin6.dat" :
             (icomp == 6) ? "input_files/iso-thin7.dat" :
             (icomp == 7) ? "input_files/iso-thick2.dat"  :
             (icomp == 8) ? "input_files/iso-bar_age.dat" :
                            "input_files/iso-NSD.dat";
     if((fp=fopen(file1,"r"))==NULL){
        printf("can't open %s\n",file1);
        exit(1);
     }
     nage = (icomp == 0) ? 3  :
            (icomp == 1) ? 18 :
            (icomp == 2) ? 21 :
            (icomp == 3) ? 21 :
            (icomp == 4) ? 41 :
            (icomp == 5) ? 41 :
            (icomp == 6) ? 61 :
            (icomp == 7) ? 2  :
            (icomp == 8) ? 27 :
                           6;
     narry = (icomp == 0) ? 465 :
             (icomp == 1) ? 581 :
             (icomp == 2) ? 1885 :
             (icomp == 3) ? 552 :
             (icomp == 4) ? 362 :
             (icomp == 5) ? 323 :
             (icomp == 6) ? 296 :
             (icomp == 7) ? 197 :
             (icomp == 8) ? 766 :
                            340;
     double **Minis, **MIcs;
     int *nMinis;
     nMinis  = (int *)calloc(nage, sizeof(double *));
     Minis  = (double **)malloc(sizeof(double *) * nage);
     MIcs   = (double **)malloc(sizeof(double *) * nage);
     for (int j=0; j<nage; j++){
       Minis[j]  = (double *)calloc(narry, sizeof(double *));
       MIcs[j]   = (double *)calloc(narry, sizeof(double *));
     }
     int iage = 0, iagest = 0, iageen = 0, nmax= 0;
     int dtau;
     dtau = (icomp == 0) ?  10 :
            (icomp == 1) ?  85 :
            (icomp == 2) ? 100 :
            (icomp == 3) ? 100 :
            (icomp == 4) ? 100 :
            (icomp == 5) ?  50 :
            (icomp == 6) ? 100 :
            (icomp == 7) ? 100 :
            (icomp == 8) ? 100 :
                           100 ;
     while (fgets(line,1000,fp) !=NULL){
        int nwords = split((char*)" ", line, words);
        if (*words[0] == '#') continue;
        if (iagest == 0) iagest = atoi(words[0]);
        if ((atoi(words[0]) - iagest) % dtau != 0) continue;
        iageen = atoi(words[0]);  
        iage = (atoi(words[0]) - iagest + 0.5)/dtau; // Gyr
        Minis[iage][nMinis[iage]] = atof(words[1]);
        MIcs[iage][nMinis[iage]] = atof(words[3]);
        nMinis[iage]++;
        if (nMinis[iage] > nmax) nmax = nMinis[iage];
     }
     fclose(fp);
     // for Ihist
     int MIst = -6, MIen = 13, Nbin = 950;
     double dI = (double) (MIen - MIst)/Nbin, pIs[960] = {};
     double gamma = 1/tSFR;  // SFR timescale, 7 Gyr
     for (int j= 0; j< iage+1;j++){
       double tau = (j*dtau + iagest)*0.01; // Gyr
       double wtSFR = (icomp == 9) ? exp(-0.5*pow((tau - mageND)/sageND, 2))   // 7 +- 1 of Gaussian
                    : (icomp == 8) ? exp(-0.5*pow((tau - mageB)/sageB, 2))   // 9 +- 1 of Gaussian
                    : (icomp <  7) ? exp(-gamma*(10-tau))  // thin
                    : 2; // Thick
       if (j == 0 || j == iage) wtSFR *= 0.5; // daikei sekibun
       if (j == 0 && icomp == 0) wtSFR *= 3; // add 0.00 Gyr - 0.05 Gyr
       double logMini = log10(Minis[j][0]); // should be ~0.09 Msun
       double PBD = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini);
       pIs[Nbin - 5] += wtSFR * PBD;
       for (int k=0; k< nMinis[j] - 1; k++){
         if (Minis[j][k+1] == 0) continue;
         double Mini1 = Minis[j][k];
         double Mini2 = (Minis[j][k+1] > 0) ? Minis[j][k+1] : 0;
         if (Mini1 > Mini2) printf ("Worning!! Mini1 > Mini2 !!!!\n");
         double MIc = 0.5 * (MIcs[j][k] + MIcs[j][k+1]);
         double logMini1 = log10(Minis[j][k]);
         double logMini2 = log10(Minis[j][k+1]);
         double P1 = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini1);
         double P2 = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini2);
         double wtM = P2 - P1;
         int pI = (MIc - MIst)/dI;
         if (pI < 0) pI = 0;
         if (pI >= Nbin) pI = Nbin - 5;
         pIs[pI] += wtSFR * wtM; 
       }
     }

     for (int pI=0;pI<=Nbin;pI++){
        double MI = MIst + pI * dI;
        // printf ("%6.2f",MI);
        if (pI>=1) {
          CumuN_MIs[icomp][pI] = 0.5*(pIs[pI] + pIs[pI-1])  + CumuN_MIs[icomp][pI-1];
        } else {
          CumuN_MIs[icomp][pI] = 0.0;
        }
     }
     for (int j=0; j<nage; j++){
       free(Minis[j]);
       free(MIcs[j]) ;
     }
     free (Minis);
     free (MIcs);
     free (nMinis);
   }

   //----------------------------------
   // Normalize cumulative distirbution
   for (int k=0; k<ncomp; k++){
   for (int j=0; j<i; j++){
      CumuN_MIs[k][j] /= CumuN_MIs[k][i-1];
      // printf ("%1d-%03d %5.2f %.6e\n",k,j,MIs[j],CumuN_MIs[k][j]);
   }}
   // printf ("#N= %d read from %s\n",i, infile);


   return i;
}
//---------------
void store_VI_MI(double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI, double *MIs, double *VIs, double ***f_VI_Is, double *logMass, double *PlogM_cum_norm){
  const char *file1;
  char line[1000];
  char *words[100];
  FILE *fp;
  int nage, narry;
  double dMI = (double) (MIen - MIst)/NbinMI;
  double dVI = (double) (VIen - VIst)/NbinVI;
  for (int icomp=0; icomp<ncomp; icomp++){
    file1 = (icomp == 0) ? "input_files/iso-thin1.dat" :
            (icomp == 1) ? "input_files/iso-thin2.dat" :
            (icomp == 2) ? "input_files/iso-thin3.dat" :
            (icomp == 3) ? "input_files/iso-thin4.dat" :
            (icomp == 4) ? "input_files/iso-thin5.dat" :
            (icomp == 5) ? "input_files/iso-thin6.dat" :
            (icomp == 6) ? "input_files/iso-thin7.dat" :
            (icomp == 7) ? "input_files/iso-thick2.dat"  :
            (icomp == 8) ? "input_files/iso-bar_age.dat" :
                           "input_files/iso-NSD.dat";
    if((fp=fopen(file1,"r"))==NULL){
       printf("can't open %s\n",file1);
       exit(1);
    }
    nage = (icomp == 0) ? 3  :
           (icomp == 1) ? 18 :
           (icomp == 2) ? 21 :
           (icomp == 3) ? 21 :
           (icomp == 4) ? 41 :
           (icomp == 5) ? 41 :
           (icomp == 6) ? 61 :
           (icomp == 7) ? 2 : 
           (icomp == 8) ? 27 :
                          6;
    narry = (icomp == 0) ? 465 :
            (icomp == 1) ? 581 :
            (icomp == 2) ? 1885 :
            (icomp == 3) ? 552 :
            (icomp == 4) ? 362 :
            (icomp == 5) ? 323 :
            (icomp == 6) ? 296 :
            (icomp == 7) ? 197 :
            (icomp == 8) ? 766 :
                           340;
    double **Minis, **VIcs, **MIcs, *ages;
    int *nMinis;
    nMinis  = (int *)calloc(nage, sizeof(int    *));
    ages    = (double *)calloc(nage, sizeof(double *));
    Minis  = (double **)malloc(sizeof(double *) * nage);
    VIcs   = (double **)malloc(sizeof(double *) * nage);
    MIcs   = (double **)malloc(sizeof(double *) * nage);
    for (int j=0; j<nage; j++){
      Minis[j]  = (double *)calloc(narry, sizeof(double *));
      VIcs[j]   = (double *)calloc(narry, sizeof(double *));
      MIcs[j]   = (double *)calloc(narry, sizeof(double *));
    }
    int iage = -1;
    int dtau;
    double age0 = 0;
    while (fgets(line,1000,fp) !=NULL){
       int nwords = split((char*)" ", line, words);
       if (*words[0] == '#') continue;
       double age = atof(words[0]);
       if (age != age0) iage++;
       age0 = age;
       Minis[iage][nMinis[iage]] = atof(words[1]);
       VIcs[iage][nMinis[iage]] = atof(words[2]) - atof(words[3]);
       MIcs[iage][nMinis[iage]] = atof(words[3]);
       if (nMinis[iage] == 0) ages[iage] = age;
       nMinis[iage]++;
    }
    fclose(fp);
    // for Ihist
    double gamma = 1/tSFR;  // SFR timescale, 7 Gyr
    double sumwt = 0;
    for (int j= 0; j< iage+1;j++){
      double tau = 0.01*ages[j]; // in Gyr
      double wtSFR = (icomp == 9) ? exp(-0.5*pow((tau - mageND)/sageND, 2))   // 7 +- 1 of Gaussian
                   : (icomp == 8) ? exp(-0.5*pow((tau - mageB)/sageB, 2))   // 9 +- 1 of Gaussian
                   : (icomp <  7) ? exp(-gamma*(10-tau))  // thin
                   : 2; // Thick
      if (j == 0 || j == iage) wtSFR *= 0.5; // daikei sekibun
      if (j == 0 && icomp == 0) wtSFR *= 3; // add 0.00 Gyr - 0.05 Gyr
      double logMini = log10(Minis[j][0]); // should be ~0.09 Msun
      double PBD = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini);
      f_VI_Is[icomp][NbinVI - 5][NbinMI - 5] += wtSFR * PBD;
      sumwt += wtSFR * PBD;
      for (int k=0; k< nMinis[j] - 1; k++){
        if (Minis[j][k+1] == 0) continue;
        double Mini1 = Minis[j][k];
        double Mini2 = (Minis[j][k+1] > 0) ? Minis[j][k+1] : 0;
        if (Mini1 > Mini2) printf ("Warning!! Mini1 > Mini2 !!!! @ tau= %.2f M1= %.9f M2= %.9f\n",tau,Mini1,Mini2);
        double logMini1 = log10(Minis[j][k]);
        double logMini2 = log10(Minis[j][k+1]);
        double P1 = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini1);
        double P2 = interp_x(nm+1, PlogM_cum_norm, logMst, dlogM, logMini2);
        double wtM = P2 - P1;
        for (int l=0; l<3; l++){
          double VI  = (l == 0) ? VIcs[j][k] 
                     : (l == 1) ? 0.5 * (VIcs[j][k] + VIcs[j][k+1])
                     : VIcs[j][k+1];
          double MI  = (l == 0) ? MIcs[j][k] 
                     : (l == 1) ? 0.5 * (MIcs[j][k] + MIcs[j][k+1])
                     : MIcs[j][k+1];
          int pMI = (MI - MIst)/dMI;
          if (pMI < 0) pMI = 0;
          if (pMI >= NbinMI) pMI = NbinMI - 1;
          int pVI = (VI - VIst)/dVI;
          if (pVI < 0) pVI = 0;
          if (pVI >= NbinVI) pVI = NbinVI - 1;
          double tmpwt = (l == 1) ? 1 : 0.5;
          f_VI_Is[icomp][pVI][pMI] += tmpwt* wtSFR * wtM; 
          sumwt += tmpwt* wtSFR * wtM;
        }
      }
    }

    for (int pVI=0;pVI<=NbinVI;pVI++){
       double VI = VIst + pVI * dVI;
       if (icomp == 0)  VIs[pVI] = VI;
       for (int pMI=0;pMI<=NbinMI;pMI++){
          double MI = MIst + pMI * dMI;
          if (icomp == 0) MIs[pMI] = MI;
          f_VI_Is[icomp][pVI][pMI] /= sumwt;
       }
    }
    for (int j=0; j<nage; j++){
      free(Minis[j]);
      free(VIcs[j]) ;
      free(MIcs[j]) ;
    }
    free (Minis);
    free (VIcs);
    free (MIcs);
    free (nMinis);
  }
}

//---- calc P(Rg|R) following Shu distribution ( Eq.(14) of Sharma et al. 2014, ApJ, 793, 51) -------
double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd){ 
  double getx2y(int n, double *x, double *y, double xin);
  if (fg <= 0) return 0;
  double  calc_faca(double Rg, double hsigU, int rd, double a0);
  double calc_SigRg(double Rg, double hsigU, int rd, double a0);
  double    calc_gc(double c);
  double Rg = R*fg;
  double vc = getx2y(nVcs, Rcs, Vcs, Rg) / (1 + 0.0374*pow(0.001*fabs(z), 1.34));
  double a0 = sigU0/vc * exp(R0/hsigU);
  double a  = sigU0/vc * exp(-(Rg - R0)/hsigU);
  double faca = calc_faca(Rg,hsigU,rd,a0);
  a *= faca;
  double c = 0.5/a/a;
  if (c <= 0.5) return 0;
  double SigRg = calc_SigRg(Rg,hsigU,rd,a0);
  double gc = calc_gc(c);
  double x = c*(2*log(fg) + 1 - fg*fg);
  double PRRg = SigRg * exp(x)/gc;
  // printf ("(calc_PRRg) R=%d, Rg=%f, vc=%f a0= sigU0(%f)/vc(%f) * exp(R0(%f)/hsigU(%f))(%.3e)= %.3e, a= %.3e, c= %.3e, SigRg= %.3e, gc= %.3e, x= %.3e, PRRg= SigRg*exp(x)/gc = %.3e * %.3e = %.3e\n",R,Rg,vc,sigU0,vc,R0,hsigU,exp(R0/hsigU),a0,a,c,SigRg,gc,x,SigRg,exp(x)/gc,PRRg);
  if (PRRg < 0) PRRg = 0;
  return PRRg;
}
/*----------------------------------------------------------------*/
double calc_gc(double c){ // Eq.(16) of Sharma et al. 2014, ApJ, 793, 51
   if (c < 0.5) return 0;
  double c2, gamma, c3, gc;
  if (c < 10){
    c2 = c - 0.5;
    gamma = tgamma(c2);
    c3 = 2 * pow(c, c2);
    gc = exp(c) * gamma/c3;
  }else{
    gc = sqrt(0.5*PI/(c - 0.913)); // approximation Eq. (14) of Schonrich & Binney 2012
  }
  return gc;
}
/*----------------------------------------------------------------*/
double calc_SigRg(double Rg, double hsigU, int rd, double a0){ // Rd**2 * Eq.(20) of Sharma et al. 2014, ApJ, 793, 51
  double k = 31.53, a = 0.6719, b = 0.2743;
  // my (c1, c2, c3, c4) = (3.740, 0.523, 0.00976, 2.29); # for flat vc
  double c1 = 3.822, c2 = 0.524, c3 = 0.00567, c4 = 2.13; // for rising vc from Table 1 of Sharma & Bland-Hawhorn (2013), ApJ, 773, 183
  double q = rd/hsigU;
  double Rgmax = c1*rd/(1+q/c2); // Eq.32 of Sharma & Bland-Hawhorn (2013), ApJ, 773, 183
  // x = Rg/3.74/Rd/(1+q/0.523); # This is form in Sharma+14, but wrong
  double x = Rg/Rgmax;  // x = Rg/Rgmax in Sharma & Bland-Hawhorn (2013), ApJ, 773, 183
  double s = k*exp(-x/b)*((x/a)*(x/a) - 1); // Eq. (21) of Sharma et al. 2014, ApJ, 793, 51
  double SigRg = 0.5*exp(-Rg/rd)/PI - c3*pow(a0,c4) * s;
  return SigRg;
}
/*----------------------------------------------------------------*/
double calc_faca(double Rg, double hsigU, int rd, double a0){ // Eq.(39) of Sharma & Bland-Hawhorn (2013), ApJ, 773, 183
  double q = rd/hsigU;
  double bunsi = 0.25*pow(a0, 2.04);
  double bumbo = pow(q, 0.49);
  double as[12] = {-0.028476,-1.4518,12.492,-21.842,19.130,-10.175,3.5214,-0.81052,0.12311,-0.011851,0.00065476,-1.5809e-05};
  double x = Rg*q/rd;
  double fpoly = as[0] + as[1]*x + as[2]*pow(x,2.) + as[3]*pow(x,3.) + as[4]*pow(x,4.) + as[5]*pow(x,5.) + as[6]*pow(x,6.) + as[7]*pow(x,7.) + as[8]*pow(x,8.) + as[9]*pow(x,9.) + as[10]*pow(x,10.) + as[11]*pow(x,11.);
  double faca = (1 - bunsi/bumbo * fpoly);
  return faca;
}


/*----------------------------------------------------------------*/
/*                   for Normalize rho or sigv                    */
/*----------------------------------------------------------------*/
double crude_integrate(double xmax, double ymax, double zmax, int nbun)  // for normalize rho_b
{
  double calc_rhoB(double xb, double yb, double zb);
  int get_p_integral(int nji, double *ls, double *ks);
  double *ls, *ks;
  int nmin, narry, nji, ncalc;
  nji  =   2;
  // nbun =  30;
  narry = (nji <= 1) ?  1 :
          (nji <= 2) ?  3 :
          (nji <= 4) ?  9 :
          (nji <= 6) ? 18 :
          (nji <= 8) ? 30 : 42;
  ls = (double *)malloc(sizeof(double *) * narry);
  ks = (double *)malloc(sizeof(double *) * narry);
  nmin = get_p_integral(nji, ls, ks);
  if (nbun < nmin) nbun = nmin;
  ncalc = nbun + 1 + 2*narry - 2*nji;  // ls[narry] includes i <= nji - 1
  double xb, xb0, yb, yb0, zb, zb0, rho, rho0, *rhosumz, *rhosumyz;
  rhosumz  = (double *)malloc(sizeof(double *) * ncalc);
  rhosumyz = (double *)malloc(sizeof(double *) * ncalc);
  double dx, dy, dz, dxtmp, dytmp, dztmp;
  // int xmax = 2200, ymax = 1400, zmax = 1200;
  dx = (double) (xmax - 0)/nbun;
  dy = (double) (ymax - 0)/nbun;
  dz = (double) (zmax - 0)/nbun;
  // printf ("dx= %.1f dy= %.1f dz= %.1f\n",dx,dy,dz);
  double totalmass = 0, massVVVbox = 0;
  int ixtmp, iytmp, iztmp;
  for (int ix = 0; ix < ncalc; ix++){
    ixtmp = ix - 2*narry + nji; // ixtmp = nji - nbun - nji
    xb = (ix>=2*narry) ? 0 + dx * ixtmp 
        : (ix % 2 == 0) ? 0 + dx * ls[ix/2] : xmax - dx * ls[ix/2];
    for (int iy = 0; iy < ncalc; iy++){
       iytmp = iy - 2*narry + nji; // iytmp = nji - nbun - nji
       yb = (iy>=2*narry) ? 0 + dy * iytmp 
          : (iy % 2 == 0) ? 0 + dy * ls[iy/2] : ymax - dy * ls[iy/2];
       rhosumz[iy] = 0;
       for(int j=0;j< narry;j++){
           dztmp = dz*ls[j];
           zb0   =    0 + dztmp;
           zb    = zmax - dztmp;
           rho0       = calc_rhoB(xb, yb,  zb0);
           rho        = calc_rhoB(xb, yb,   zb);
           rhosumz[iy] += (rho0 + rho)*ks[j];
       }
       for(int j=nji;j<=nbun-nji;j++){
           zb  =    0 + dz*j;
           rho = calc_rhoB(xb, yb, zb);
           rhosumz[iy] += rho;
       }
       rhosumz[iy] *= dz;
       // printf ("iy=%d iytmp=%d yb= %f ndy= %f ls[%d]= %f rhosumz= %f\n",iy,iytmp,yb,yb/dy,iy/2,ls[iy/2],rhosumz[iy]);
    }
    rhosumyz[ix] = 0;
    for(int j=0;j< narry;j++){
        rho0  = rhosumz[2*j];
        rho   = rhosumz[2*j+1];
        // printf ("iy=%d rhosumz= %f\n",2*j  ,rhosumz[2*j]);
        // printf ("iy=%d rhosumz= %f\n",2*j+1,rhosumz[2*j+1]);
        rhosumyz[ix] += (rho0 + rho)*ks[j];
    }
    for(int j=2*narry;j<ncalc;j++){
        rhosumyz[ix] += rhosumz[j];
    }
    rhosumyz[ix] *= dy;
    // printf ("ix=%d ixtmp=%d ndx= %f ls[%d]= %f rhosumyz= %f\n",ix,ixtmp,xb/dx,ix/2,ls[ix/2],rhosumyz[ix]);
  }
  for(int j=0;j< narry;j++){
      rho0  = rhosumyz[2*j];
      rho   = rhosumyz[2*j+1];
      // printf ("ix=%d rhosumyz= %f\n",2*j  ,rhosumyz[2*j]);
      // printf ("ix=%d rhosumyz= %f\n",2*j+1,rhosumyz[2*j+1]);
      totalmass += (rho0 + rho)*ks[j];
  }
  for(int j=2*narry;j<ncalc;j++){
      totalmass += rhosumyz[j];
  }
  totalmass *= dx*8;
  massVVVbox = totalmass;
  free (ls);
  free (ks);
  free (rhosumz);
  free (rhosumyz);
  return totalmass;
}
//---------------
void calc_sigvb(double xb, double yb, double zb, double *sigvbs)
{
  double xn, yn, zn, Rs, rs, facsig, facsigz = 0;
  xn = fabs(xb/x0_vb), yn = fabs(yb/y0_vb), zn = fabs(zb/z0_vb);
  Rs = pow((pow(xn, C1_vb) + pow(yn, C1_vb)), 1/C1_vb);
  rs = pow(pow(Rs, C2_vb) + pow(zn, C2_vb), 1/C2_vb);
  if (rs==0 && model_vb == 8) rs = 0.0001; // to avoid infty
  facsig = (model_vb == 5) ? exp(-rs)  // exponential
         : (model_vb == 6) ? exp(-0.5*rs*rs) // Gaussian
         : (model_vb == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
         : (model_vb == 4) ? exp(-pow(rs,C3_vb))  
         : 0;
  if (model_vbz >= 4){
    xn = fabs(xb/x0_vbz), yn = fabs(yb/y0_vbz), zn = fabs(zb/z0_vbz);
    Rs = pow((pow(xn, C1_vbz) + pow(yn, C1_vbz)), 1/C1_vbz);
    rs = pow(pow(Rs, C2_vbz) + pow(zn, C2_vbz), 1/C2_vbz);
    if (rs==0 && model_vbz == 8) rs = 0.0001; // to avoid infty
    facsigz = (model_vbz == 5) ? exp(-rs)  // exponential
            : (model_vbz == 6) ? exp(-0.5*rs*rs) // Gaussian
            : (model_vbz == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
            : (model_vbz == 4) ? exp(-pow(rs,C3_vbz))  
            : 0;
  }else{
    facsigz = facsig;
  }
  sigvbs[0] = sigx_vb * facsig + sigx_vb0;
  sigvbs[1] = sigy_vb * facsig + sigy_vb0;
  sigvbs[2] = sigz_vb * facsigz + sigz_vb0;
}
/*----------------------------------------------------------------*/
/*                      for general use                           */
/*----------------------------------------------------------------*/
//-----get c= a x b------------------
void cross(double *c, double *a, double *b){
   c[0] = a[1]*b[2] - b[1]*a[2];
   c[1] = a[2]*b[0] - b[2]*a[0];
   c[2] = a[0]*b[1] - b[0]*a[1];
}
//-----get c= a . b------------------
double dot(double *a, double *b){
   return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
//---- normalize a vector ----
void norm_vec(double *a){
  double norm = dot(a, a);
  for (int i=0;i<3;i++){
    a[i] = a[i]/sqrt(norm);
  }
}
//--------
double calc_rho_n(double D, int idata, double *rho_n){  // return rho, n, and wtDBs 
  double m_idisk = 0, f_disk = 0;
	// Calc rho & n
  void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb);  // return rho for each component 
  double *rhos, xyz[3] = {}, xyb[2] = {};
  rhos        = (double *)calloc(ncomp, sizeof(double *));
	calc_rho_each(D, idata, rhos, xyz, xyb);
  rho_n[0] = rho_n[1] = 0;
  if (DISK > 0) {
    for (int i=0; i<8; i++){
			rho_n[0] += rho0d[i] * rhos[i];
			rho_n[1] += n0MSd[i] * rhos[i];
      m_idisk += i * n0MSd[i] * rhos[i];
		}
    m_idisk /= rho_n[1]; // mean of idisk
    m_idisk *= 10; // 10 x mean of idisk
	}
	double nD = rho_n[1]; // number of disk stars
  // if (LB   > 0){ rho_n[0] += rho0LB   * rhos[9];
	//                rho_n[1] += n0MSLB   * rhos[9];}
  // if (STB  > 0){ rho_n[0] += rho0STB  * rhos[10];
 	//                rho_n[1] += n0MSSTB  * rhos[10];}
  // Bar
  rho_n[0] += rho0b * rhos[8] + rho0ND * rhos[9];
  rho_n[1] += n0MSb * rhos[8] + n0MSND * rhos[9];
  f_disk = nD/rho_n[1];
  free(rhos);
  return (int) m_idisk + 0.1*f_disk;
}
//---------------
void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb){  // return rho for each component 
  void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);
  double calc_rhoB(double xb, double yb, double zb);
  double x, y, z, R, xb, yb, zb, xn, yn, zn, rs, zdtmp, rhotmp;
  double lD, bD;
  lD = lDs[idata];
  bD = bDs[idata];
  Dlb2xyz(D, lD, bD, R0, xyz);
  x = xyz[0], y = xyz[1], z = xyz[2];
  R = sqrt(x*x + y*y);
  // i = 0-6: thin disk, i=7: thick disk, i=8: bulge, i=9: long bar, i = 10: super thin bar
  for (int i = 0; i<ncomp; i++){rhos[i] = 0;} // shokika
  // Disk
  int idisk, itmp, ist;
  if (DISK > 0){
    double ftmp = (hDISK == 0) ? 0.005 : (hDISK == 1) ? 0.01 : 0;
    // ist = (fabs(z) < 400) ? ftmp*fabs(z)  : (fabs(z) <= 1200) ?  4 : 7;
    ist = 0;
    for (idisk = ist; idisk < 8; idisk++){ // ignore disk0 - disk3
      zdtmp =(hDISK == 0) ? zd[idisk] :
              (R > 4500)  ? zd[idisk] + (R-R0)*(zd[idisk] - zd45[idisk])/(R0 - 4500) 
                          : zd45[idisk];
      rhotmp  = (idisk < 7) ? 4.0/(exp(2*z/zdtmp)+exp(-2*z/zdtmp)+2)
                            : exp(-fabs(z)/zd[idisk]);
      itmp = (idisk == 0) ? 0 : (idisk <  7) ? 1 : 2;
      rhotmp *= zd[idisk]/zdtmp;  // zd/zdtmp is to keep Sigma(R) as exponential 
      if (DISK == 1) rhotmp = rhotmp * exp(-R/Rd[itmp] - pow(((double)Rh/R),nh));
      if (DISK == 2) rhotmp = (R > Rdbreak) ? rhotmp * exp(-R/Rd[itmp])
                                            : rhotmp * exp(- (double)Rdbreak/Rd[itmp]); // const. in R < 5300
      if (DISK == 3) rhotmp = rhotmp * exp(-R/Rd[itmp]);
      rhos[idisk]  = rhotmp/y0d[itmp]; // Number density of BD + MS
    }
  }
  // Bar
  xb =  x * costheta + y * sintheta;
  yb = -x * sintheta + y * costheta;
  zb =  z;                          
  rhos[8] = calc_rhoB(xb,yb,zb);
  // ND 
  if (ND > 0){
    if (ND == 3){
      if (R <= RenND - 30 && fabs(z) <= zenND - 20){
        rhos[9] = pow(10.0, interp_xy(nzND, nRND, logrhoNDs, zstND, RstND, dzND, dRND, fabs(z), R));
      }else{
        rhos[9] = 0;
      }
    }else{
      // See Eq. (28) of Portail et al. 2017
      xn = fabs(xb/x0ND), yn = fabs(yb/y0ND), zn = fabs(zb/z0ND);
      rs = pow((pow(xn, C1ND) + pow(yn, C1ND)), 1/C1ND) + zn;
      rhos[9] = exp(-rs);  
    }
  }
  xyb[0] = xb;
  xyb[1] = yb;
}
//---------------
double calc_rhoB(double xb, double yb, double zb)
{
  double xn, yn, zn, R, Rs, rs, rho, rho2, rhoX;
   R = sqrt(xb*xb + yb*yb);
  // 1st  Bar
  if (model >= 4 && model <= 8){
    xn = fabs(xb/x0_1), yn = fabs(yb/y0_1), zn = fabs(zb/z0_1);
    Rs = pow((pow(xn, C1) + pow(yn, C1)), 1/C1);
    rs = pow(pow(Rs, C2)     + pow(zn, C2), 1/C2);
    if (rs==0 && model == 8) rs = 0.0001; // to avoid infty
    rho = (model == 5) ? exp(-rs)  // exponential for 4 or 5
        : (model == 6) ? exp(-0.5*rs*rs) // Gaussian
        : (model == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
        : (model == 4) ? exp(-pow(rs, C3))
        : 0;
  }
  if (R  >= Rc) rho *= exp(-0.5*(R-Rc)*(R-Rc)/srob/srob);
  if (fabs(zb) >= zb_c) rho *= exp(-0.5*(fabs(zb)-zb_c)*(fabs(zb)-zb_c)/200.0/200.0);

  // X-shape
  if (addX >= 5){
    xn = fabs((xb-b_zX*zb)/x0_X), yn = fabs((yb-b_zY*zb)/y0_X), zn = fabs(zb/z0_X);
    rs = pow(pow((pow(xn, C1_X) + pow(yn, C1_X)), C2_X/C1_X) + pow(zn, C2_X), 1/C2_X);
    rhoX  = (addX == 5) ? exp(-rs)  // exponential
          : (addX == 6) ? exp(-0.5*rs*rs) // Gaussian
          : (addX == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
           : 0;
    xn = fabs((xb+b_zX*zb)/x0_X), yn = fabs((yb-b_zY*zb)/y0_X);
    rs = pow(pow((pow(xn, C1_X) + pow(yn, C1_X)), C2_X/C1_X) + pow(zn, C2_X), 1/C2_X);
    rhoX += (addX == 5) ? exp(-rs)  // exponential
          : (addX == 6) ? exp(-0.5*rs*rs) // Gaussian
          : (addX == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
           : 0;
    if (b_zY > 0.0){
      xn = fabs((xb-b_zX*zb)/x0_X), yn = fabs((yb+b_zY*zb)/y0_X);
      rs = pow(pow((pow(xn, C1_X) + pow(yn, C1_X)), C2_X/C1_X) + pow(zn, C2_X), 1/C2_X);
      rhoX += (addX == 5) ? exp(-rs)  // exponential
            : (addX == 6) ? exp(-0.5*rs*rs) // Gaussian
            : (addX == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
             : 0;
      xn = fabs((xb+b_zX*zb)/x0_X), yn = fabs((yb+b_zY*zb)/y0_X);
      rs = pow(pow((pow(xn, C1_X) + pow(yn, C1_X)), C2_X/C1_X) + pow(zn, C2_X), 1/C2_X);
      rhoX += (addX == 5) ? exp(-rs)  // exponential
            : (addX == 6) ? exp(-0.5*rs*rs) // Gaussian
            : (addX == 7) ? pow( 2.0/(exp(rs)+exp(-rs)), 2) // sech2
             : 0;
    }
    rhoX *= fX;
    if (R >= Rc_X) rhoX *= exp(-0.5*(R-Rc_X)*(R-Rc_X)/srob/srob);
  }
  if (addX >=5) rho += rhoX;
  return rho;
}
//---------------
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz)
/*------------------------------------------------------------*/
/*  Give (x,y,z) for a given D, lD, bD  */
/*------------------------------------------------------------*/
{
  double cosbsun = cos(zsun/Rsun), sinbsun = sin(zsun/Rsun);
  double cosb = cos(bD/180.0*PI), sinb = sin(bD/180.0*PI), 
         cosl = cos(lD/180.0*PI), sinl = sin(lD/180.0*PI);
  double xtmp = Rsun - D * cosb * cosl;
  double ytmp = D * cosb * sinl;       
  double ztmp = D * sinb;              
  // xyz[0] = -ztmp * sinbsun + xtmp * cosbsun;
  xyz[0] =  xtmp - xyzSgrA[0]; 
  xyz[1] =  ytmp - xyzSgrA[1];                            
  xyz[2] =  ztmp * cosbsun + xtmp * sinbsun - xyzSgrA[2]; 
}


/*----------------------------------------------------------------*/
/*                      for general use                           */
/*----------------------------------------------------------------*/
//---- get parameters for trapezoidal integral -------
// Values from http://midarekazu.g2.xrea.com/Newton-Cotes.html
int get_p_integral(int nji, double *ls, double *ks)
{
   int nmin, i = 0, j=0, k=0;
   nji = (nji <= 1) ? 1 :
         (nji <= 2) ? 2 :
         (nji <= 4) ? 4 :
         (nji <= 6) ? 6 :
         (nji <= 8) ? 8 : 10;
   if (nji == 1){ // Trapezoid
      ls[i++]= 0.0;
      ks[j++]= 0.5;
      nmin = 1;
   }
   if (nji == 2){ // simpson
      ls[i++]= 0.0,    ls[i++]= 0.5,    ls[i++]=  1.0;
      ks[j++]= 3.0/12, ks[j++]= 4.0/12, ks[j++]= 11.0/12;
      nmin = 3;
   }
   if (nji == 4){ // boole
      ls[i++]=0.0,
      ls[i++]=1.0/4, ls[i++]=1.0/2, ls[i++]=3.0/4, ls[i++]=1.0, ls[i++]=3.0/2, ls[i++]=2.0, 
      ls[i++]=9.0/4, ls[i++]=3.0;
      ks[j++]= 70./360,  
      ks[j++]= 32./360, ks[j++]= 76./360, ks[j++]= 128./360, ks[j++]=187./360, ks[j++]= 100./360, 
      ks[j++]=218./360, ks[j++]= 96./360, ks[j++]= 353./360; 
      nmin = 7;
   }
   if (nji == 6){
      ls[i++]=0.0,
      ls[i++]=1.0/6, ls[i++]=1.0/3, ls[i++]=1.0/2, ls[i++]=2.0/3, ls[i++]=5.0/6, ls[i++]=1.0, 
      ls[i++]=4.0/3, ls[i++]=3.0/2, ls[i++]=5.0/3, ls[i++]=2.0,   ls[i++]=5.0/2, ls[i++]=8.0/3, 
      ls[i++]=3.0,   ls[i++]=10.0/3,ls[i++]=4.0,   ls[i++]=25.0/6,ls[i++]=5.0;
      ks[j++]= 861./5040,  
      ks[j++]= 216./5040,  ks[j++]= 459./5040,  ks[j++]= 920./5040,  ks[j++]= 945./5040,  
      ks[j++]=1296./5040,  ks[j++]=2208./5040,  ks[j++]= 162./5040,  ks[j++]= 816./5040,   
      ks[j++]= 567./5040,  ks[j++]=2955./5040,  ks[j++]=2008./5040,  ks[j++]= 108./5040,   
      ks[j++]=3459./5040,  ks[j++]= 999./5040,  ks[j++]=3662./5040,  ks[j++]=1080./5040,  
      ks[j++]=4999./5040;
      nmin = 11;
   }
   if (nji == 8){
      ls[i++]=0.0,
      ls[i++]=1.0/8,  ls[i++]=1.0/4, ls[i++]=3.0/8, ls[i++]=1.0/2, ls[i++]=5.0/8,  ls[i++]=3.0/4, 
      ls[i++]=7.0/8,  ls[i++]=1.0,   ls[i++]=9.0/8, ls[i++]=5.0/4, ls[i++]=3.0/2,  ls[i++]=7.0/4, 
      ls[i++]=15.0/8, ls[i++]=2.0,   ls[i++]=9.0/4, ls[i++]=5.0/2, ls[i++]=21.0/8, ls[i++]=3.0,  
      ls[i++]=25.0/8, ls[i++]=7.0/2, ls[i++]=15.0/4,ls[i++]=4.0,   ls[i++]=35.0/8, ls[i++]=9.0/2, 
      ls[i++]=5.0,    ls[i++]=21.0/4,ls[i++]=6.0,   ls[i++]=49.0/8,ls[i++]=7.0;
      ks[j++]= 35604./226800,  
      ks[j++]=  5888./226800, ks[j++]= 10848./226800, ks[j++]= 28160./226800, 
      ks[j++]= 17156./226800, ks[j++]= 39936./226800, ks[j++]= 52608./226800, 
      ks[j++]= 47104./226800, ks[j++]= 43213./226800, ks[j++]= 31488./226800, 
      ks[j++]= 16352./226800, ks[j++]= 20940./226800, ks[j++]=  5280./226800, 
      ks[j++]= 83968./226800, ks[j++]= 31410./226800, ks[j++]= 60192./226800, 
      ks[j++]= 19284./226800, ks[j++]= 91136./226800, ks[j++]=103575./226800, 
      ks[j++]= 52480./226800, ks[j++]= -8228./226800, ks[j++]= 58336./226800, 
      ks[j++]= 99196./226800, ks[j++]=102912./226800, ks[j++]= -5568./226800, 
      ks[j++]=184153./226800, ks[j++]= 28832./226800, ks[j++]=177718./226800, 
      ks[j++]= 41216./226800, ks[j++]=225811./226800;
      nmin = 15;
   }
   if (nji == 10){
      ls[i++]=0.0,
      ls[i++]=1.0/10, ls[i++]=1.0/5,  ls[i++]=3.0/10, ls[i++]=2.0/5,  ls[i++]=1.0/2,  
      ls[i++]=3.0/5,  ls[i++]=7.0/10, ls[i++]=4.0/5,  ls[i++]=9.0/10, ls[i++]=1.0,   
      ls[i++]=6.0/5,  ls[i++]=7.0/5,  ls[i++]=3.0/2,  ls[i++]=8.0/5,  ls[i++]=9.0/5, 
      ls[i++]=2.0,    ls[i++]=21.0/10,ls[i++]=12.0/5, ls[i++]=5.0/2,  ls[i++]=27.0/10, 
      ls[i++]=14.0/5, ls[i++]=3.0,    ls[i++]=16.0/5, ls[i++]=7.0/2,  ls[i++]=18.0/5, 
      ls[i++]=4.0,    ls[i++]=21.0/5, ls[i++]=9.0/2,  ls[i++]=24.0/5, ls[i++]=49.0/10,
      ls[i++]=5.0,    ls[i++]=27.0/5, ls[i++]=28.0/5, ls[i++]=6.0,    ls[i++]=63.0/10, 
      ls[i++]=32.0/5, ls[i++]=7.0,    ls[i++]=36.0/5, ls[i++]=8.0,    ls[i++]=81.0/10, 
      ls[i++]=9.0;
      ks[j++]=  883685./5987520,
      ks[j++]=  106300./5987520, ks[j++]=  164075./5987520, ks[j++]=  591300./5987520, 
      ks[j++]=   67600./5987520, ks[j++]=  958868./5987520, ks[j++]=  776475./5987520, 
      ks[j++]= 1016500./5987520, ks[j++]=   86675./5987520, ks[j++]= 1880200./5987520, 
      ks[j++]= 1851848./5987520, ks[j++]= -504300./5987520, ks[j++]=  205125./5987520, 
      ks[j++]= 2644104./5987520, ks[j++]=-1527450./5987520, ks[j++]=  628625./5987520, 
      ks[j++]= 1177276./5987520, ks[j++]= 2724000./5987520, ks[j++]= -571875./5987520, 
      ks[j++]= 2136840./5987520, ks[j++]= 2770500./5987520, ks[j++]= -734250./5987520, 
      ks[j++]= 4772079./5987520, ks[j++]=-2278500./5987520, ks[j++]= 4353576./5987520, 
      ks[j++]=-3483050./5987520, ks[j++]= 4097507./5987520, ks[j++]= -189450./5987520, 
      ks[j++]= 4377812./5987520, ks[j++]=-2375550./5987520, ks[j++]= 1906800./5987520, 
      ks[j++]= 5210935./5987520, ks[j++]=-1707150./5987520, ks[j++]= 1839525./5987520, 
      ks[j++]= 2621502./5987520, ks[j++]= 3195700./5987520, ks[j++]= -388200./5987520, 
      ks[j++]= 5361569./5987520, ks[j++]=  413675./5987520, ks[j++]= 4892386./5987520, 
      ks[j++]=  956700./5987520, ks[j++]=5971453./5987520;
      nmin = 19;
   }
   return nmin;
}
/*----------------------------------------------------------------*/
//---- getx2y for linear interpolation
double getx2y(int n, double *x, double *y, double xin)
{
   int i;
   double xmin,xmax;
   if (x[0] < x[n-1]){xmin=x[0];xmax=x[n-1];}else{xmin=x[n-1];xmax=x[0];}
   //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);

   if (xmin > xin) return 0;
   if (xmax < xin) return 0;
   //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);

   for(i=0;i<n;i++){
      //printf("i= %d x= %f logM= %f\n",i, x[i],logM);
      if (i == 0) continue;
      if (x[i] <= xin && x[i-1] >=xin || x[i] >=xin && x[i-1] <= xin){
         double yreq = (y[i]-y[i-1])/(x[i]-x[i-1])*(xin -x[i-1]) + y[i-1];
         return yreq;
      }
   }
   return 0;
}
//---- getx2y_ist for linear interpolation
double getx2y_ist(int n, double *x, double *y, double xin, int *ist)
{
   int i;
   /* The followings are commented cuz Mag vs Mini  */
   // double xmin,xmax;
   // if (x[0] < x[n-1]){xmin=x[0];xmax=x[n-1];}else{xmin=x[n-1];xmax=x[0];}
   // //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);

   // if (xmin > xin) return 0;
   // if (xmax < xin) return 0;
   // //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);

   for(i=*ist;i<n;i++){
      //printf("i= %d x= %f logM= %f\n",i, x[i],logM);
      if (i == 0) continue;
      if (x[i] <= xin && x[i-1] >=xin || x[i] >=xin && x[i-1] <= xin){
         double yreq = (y[i]-y[i-1])/(x[i]-x[i-1])*(xin -x[i-1]) + y[i-1];
         *ist = i;
         return yreq;
      }
   }
   return 0;
}
//---- getx2y_khi for linear interpolation
double getx2y_khi(int n, double *x, double *y, double xin, int *khi)
{
   int i;
   double xmin,xmax;
   if (x[0] < x[n-1]){xmin=x[0];xmax=x[n-1];}else{xmin=x[n-1];xmax=x[0];}
   //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);

   if (xmin > xin) return 0;
   if (xmax < xin) return 0;
   //printf("n=%d %f %f %f\n",n, xmin,xmax,logM);
   int klo;
   if (*khi > 0){
     klo = *khi - 1;
   }else{
     klo = 0;
     *khi = n-1;
     while (*khi-klo > 1) {
       int k=(*khi+klo) >> 1;
       if (x[k] > xin) *khi=k;
       else klo=k;
     }
   }
   double h = x[*khi] - x[klo];
   if (h == 0.0) return 0;
   double a = (x[*khi]-xin)/h;
   double b = (xin-x[klo])/h;
   double yreq = a*y[klo] + b*y[*khi];
   return yreq;
}
//---------------
double interp_x(int n, double *F, double xst, double dx, double xreq) // just for this code
{
  int    ix   = (xreq - xst)/dx;
  double xres = (xreq - xst)/dx - ix;
  if (ix < 0 || ix > n - 1) return 0;
  if (ix+1 > n - 1) return F[ix];
  double F1 = F[ix];
  double F2 = F[ix+1];
  return F1 * (1 - xres) + F2 * xres;
}
//---------------
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq) // just for this code
{
  int    ix   = (xreq - xst)/dx;
  double xres = (xreq - xst)/dx - ix;
  int    iy   = (yreq - yst)/dy;
  double yres = (yreq - yst)/dy - iy;
  if (ix < 0 || ix > nx - 1 || iy < 0 || iy > ny -1) return 0;
  if (ix+1 > nx - 1 && iy+1 > ny - 1) return F[ix][iy]; // return edge value
  if (ix+1 > nx - 1) return F[ix][iy] * (1 - yres) + F[ix][iy+1] * yres; // only interpolate y
  if (iy+1 > ny - 1) return F[ix][iy] * (1 - xres) + F[ix+1][iy] * xres; // only interpolate x
  double a1 = (1 - xres) * (1 - yres), F1 = F[ix][iy]    ;
  double a2 =      xres  * (1 - yres), F2 = F[ix+1][iy]  ;
  double a3 = (1 - xres) *      yres , F3 = F[ix][iy+1]  ;
  double a4 =      xres  *      yres , F4 = F[ix+1][iy+1];
  return a1 * F1 + a2 * F2 + a3 * F3 + a4 * F4;
}
//---------------
void interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq) // just for this code
/* Return coefficients as[0], as[1], as[2], as[3] of bilinear interpolation.
 * Interpolated value is given by
 *   as[0]*F[ix][iy] + as[1]*F[ix+1][iy] + as[2]*F[ix][iy+1] + as[3]*F[ix+1][iy+1]  */
{
  int    ix   = (xreq - xst)/dx;
  double xres = (xreq - xst)/dx - ix;
  int    iy   = (yreq - yst)/dy;
  double yres = (yreq - yst)/dy - iy;
  if (ix < 0 || ix > nx - 1 || iy < 0 || iy > ny -1){ // Out of range
    as[0] = as[1] = as[2] = as[3] = 0;
  }else if (ix+1 > nx - 1 && iy+1 > ny - 1){ // return edge value
    as[1] = as[2] = as[3] = 0;
    as[0]= 1; 
  }else if (ix+1 > nx - 1){ // only interpolate y
    as[0] = (1 - yres);
    as[2] = yres;
    as[1] = as[3] = 0;
  }else if (iy+1 > ny - 1){ // only interpolate x
    as[0] = (1 - xres);
    as[1] = xres;
    as[2] = as[3] = 0;
  }else{ // when either ix or iy is not at the edge
    as[0] = (1 - xres) * (1 - yres);
    as[1] =      xres  * (1 - yres);
    as[2] = (1 - xres) *      yres ;
    as[3] =      xres  *      yres ;
  }
}

