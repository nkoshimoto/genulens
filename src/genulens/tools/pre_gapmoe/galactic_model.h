#ifndef GALACTIC_MODEL_H
#define GALACTIC_MODEL_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// ---- Physical / unit constants ----
#define PI       3.1415926535897932385
#define KS2MY    210.949526569698696   // km/s/pc -> mas/yr
#define GC       4.30091e-03           // pc Msun^-1 (km/s)^2
#define zsun     25.0                  // pc
#define srob     500.0                 // pc (Gaussian cutoff width for bar)
#define MNSMIN   1.2                   // Msun
#define MNSMAX   2.1                   // Msun

// ---- GSL RNG ----
extern const gsl_rng_type *T_rng;
extern gsl_rng *r_rng;

double ran1();
double gasdev();

// ---- Component counts ----
// 0-6: thin disk (7 age bins), 7: thick disk, 8: bar, 9: NSD, 10: stellar halo
extern int ncomp;   // = 11

// ---- Global model flags ----
extern int DISK;    // 0: no disk, 1: w/ hole, 2: Portail+17, 3: simple exp
extern int hDISK;   // 0: const scale height, 1: linear in R
extern int addX;    // 0: no X-shape, >=5: add X-bar model
extern int model;   // bulge shape: 4-8
extern int ND;      // NSD: 0: none, 1: P17, 2: Sormani-like, 3: moments
extern int SH;      // stellar halo: 0: none, 1: Robin+03

// ---- Bar geometry ----
extern double R0;
extern double thetaD;
extern double x0_1, y0_1, z0_1;
extern double C1, C2, C3;
extern double Rc, frho0b, zb_c;
extern double costheta, sintheta;

// ---- X-shape bar ----
extern double x0_X, y0_X, z0_X;
extern double C1_X, C2_X;
extern double b_zX, b_zY, fX, Rc_X;

// ---- Disk scale lengths/heights ----
extern int    Rd[3];     // scale lengths (pc): thin1, thin2-7, thick
extern int    Rh, Rdbreak, nh;
extern double zd[8];     // scale heights (pc)
extern double zd45[8];   // scale heights at R=4500 pc (for hDISK=1)
extern double y0d[3];    // disk density normalization at solar position

// ---- Disk SFR ----
extern double tSFR;      // Gyr
extern double medtauds[8];  // median ages of disk populations (Gyr)

// ---- Disk kinematics ----
extern double hsigUt, hsigWt, hsigUT, hsigWT;
extern double betaU, betaW, sigU10d, sigW10d, sigU0td, sigW0td;

// ---- Bar kinematics ----
extern int    model_vb, model_vbz;
extern double Omega_p;
extern double x0_vb, y0_vb, z0_vb, C1_vb, C2_vb, C3_vb;
extern double sigx_vb, sigy_vb, sigz_vb;
extern double sigx_vb0, sigy_vb0, sigz_vb0;
extern double vx_str, y0_str;
extern double x0_vbz, y0_vbz, z0_vbz, C1_vbz, C2_vbz, C3_vbz;

// ---- Sun kinematics ----
extern double vxsun, Vsun, vzsun, vysun;

// ---- Density normalizations ----
extern double rhot0;        // local thin disk mass density (Msun/pc^3)
extern double rho0d[8];     // mass density at solar pos (Msun/pc^3)
extern double n0d[8];       // number density (MS+WD) at solar pos
extern double n0MSd[8];     // MS number density at solar pos
extern double n0RGd[8];     // RG number density at solar pos

extern double rho0b, n0MSb, n0RGb, n0b;  // bar normalizations
extern double fb_MS, m2nb_MS, m2nb_WD, nMS2nRGb;

// ---- NSD ----
extern int    x0ND, y0ND, z0ND;
extern double C1ND, rho0ND, n0MSND, n0RGND, n0ND;
extern double fND_MS, m2nND_MS, m2nND_WD, nMS2nRGND;
// NSD moment tables (ND==3)
extern double **logrhoNDs, **vphiNDs, ***logsigvNDs, **corRzNDs;
extern double zstND, zenND, dzND;
extern double RstND, RenND, dRND;
extern int    nzND, nRND;

// ---- Stellar halo ----
extern double rho0SHMS, epsSH, alphaSH, acSH2;
extern double rho0SH, n0MSSH, n0RGSH, n0SH;
extern double sigU_SH, sigV_SH, sigW_SH;

// ---- Age/death-mass tables (set by store_IMF_nBs from Minidie.dat) ----
extern int    agesD[250], agesB[50], agesND[10];
extern double MinidieD[250], MinidieB[50], MinidieND[10];
extern int    nageD, nageB, nageND;
extern double MiniWDmax;
extern double mageB, sageB, mageND, sageND;

// ---- IMF ----
extern int    nm;
extern double logMst, dlogM;

// ---- Shu DF kinematic tables (Phase 3, allocated only when needed) ----
extern double ****fgsShu, ****PRRgShus, ****cumu_PRRgs;
extern int   ***n_fgsShu, ****kptiles;
extern int    zstShu, zenShu, dzShu;
extern int    RstShu, RenShu, dRShu;

// ---- Circular velocity table (for kinematics) ----
extern int    nVcs;
extern double Rcs[60], Vcs[60];

// ---- Line-of-sight arrays ----
extern double *lDs, *bDs;  // (l, b) per data point (we use idata=0)

// ---- SgrA* offset ----
extern double xyzSgrA[3];

// ---- B14/Koshimoto+21 model flags ----
extern int B14disk, B14vbar;

// ---- IMF table (retained after init_galactic_model for tool use) ----
extern double *g_logMass;
extern double *g_PlogM;
extern double *g_PlogM_cum_norm;

// =========================================================
// Function prototypes
// =========================================================

// --- Utility ---
double getx2y(int n, double *x, double *y, double xin);
double getx2y_ist(int n, double *x, double *y, double xin, int *ist);
double getx2y_khi(int n, double *x, double *y, double xin, int *khi);
double interp_x(int n, double *F, double xst, double dx, double xreq);
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq);
void   interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq);

// --- Linear algebra ---
void   cross(double *c, double *a, double *b);
double dot(double *a, double *b);
void   norm_vec(double *a);

// --- Coordinates ---
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);
void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);

// --- Density ---
double calc_rhoB(double xb, double yb, double zb);
void   calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb);

// --- Integration ---
int    get_p_integral(int nji, double *ls, double *ks);
double crude_integrate(double xmax, double ymax, double zmax, int nbun);

// --- Data loaders ---
void store_NSDmoments(char *infile);
void store_IMF_nBs(int B, double *logMass, double *PlogM, double *PlogM_cum_norm,
                   int *imptiles,
                   double M0, double M1, double M2, double M3, double Ml, double Mu,
                   double alpha1, double alpha2, double alpha3, double alpha4, double alpha0);
void Mini2Mrem(double *pout, double Mini, int mean);

// --- Kinematics (Phase 3 -- galactic_kinematics.cpp) ---
void store_cumuP_Shu(char *infile);
void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD);

// --- Model initializer ---
// need_kinematics: 1 -> also set up Shu DF velocity tables (needed for calc_murel_dist)
void init_galactic_model(int argc, char **argv, int need_kinematics);

#endif // GALACTIC_MODEL_H
