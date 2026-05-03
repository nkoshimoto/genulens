#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "genulens/cli/option.h"
#include "genulens/io/input_data.hpp"
#include "genulens/math/interpolation.hpp"
#include "genulens/math/quadrature.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/simulation/observation_likelihood.hpp"
#include "genulens/simulation/run_context.hpp"

#define fopen(path, mode) genulens::open_input_file((path), (mode))

namespace genulens {
extern RunContext *active_state;
double ran1();
double gasdev();

struct PopulationRuntime {
  double *log_mass = nullptr;
  double *mass_probability = nullptr;
  double *mass_cumulative = nullptr;
  int *mass_percentiles = nullptr;
  double *empirical_masses = nullptr;
  double **empirical_magnitudes = nullptr;
  int empirical_count = 0;
  int empirical_columns = 5;
  bool has_luminosity_function = false;
  bool has_color_magnitude_function = false;

  void initialize_mass_function(const InitialMassFunctionOptions &options);
  void initialize_luminosity_functions(
      double source_i_min,
      double source_i_max,
      double source_vi_min,
      double source_vi_max,
      double ai_rc,
      double evi_rc);
  void read_empirical_mass_luminosity();
  void release_luminosity_functions();
  void release_all();
};

struct KinematicRuntimeTables {
  int n_z = 0;
  int n_R = 0;
  int n_disk = 8;
  int n_fg = 100;

  void initialize_shu_distribution();
  void release_all();
};

void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
void calc_opticaldepth(double *tauall, double *Nsall, int idata, int Dsmax21, double AI0, double hscale, double Isst, double Isen);
void store_NSDmoments(char *infile);
double like_obs(double mod, double obs, double err, double fe, int det, int UNIFORM);
void store_IMF_nBs(int B, double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles, double M0, double M1, double M2, double M3, double Ml, double Mu, double alpha1, double alpha2, double alpha3, double alpha4, double alpha0);
void Mini2Mrem(double *pout, double Mini, int mean);
double fLF_detect(double extI, double Imin, double Imax, int idisk);
double fIVI_detect(double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk);
void store_cumuP_Shu(char *infile);
void get_PRRGmax2(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
void calc_dpdfg(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD);
void getaproj(double *pout, double M1, double M2, int coeff);
double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
int read_MLemp(char *infile, double *M_emps, double **Mag_emps);
int make_LFs(double *MIs_arg, double **CumuN_MIs_arg, double *logMass, double *PlogM_cum_norm);
void store_VI_MI(double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI, double *MIs_arg, double *VIs_arg, double ***f_VI_Is_arg, double *logMass, double *PlogM_cum_norm);
double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd);
double calc_gc(double c);
double calc_SigRg(double Rg, double hsigU, int rd, double a0);
double calc_faca(double Rg, double hsigU, int rd, double a0);
double crude_integrate(double xmax, double ymax, double zmax, int nbun);
void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
void cross(double *c, double *a, double *b);
double dot(double *a, double *b);
void norm_vec(double *a);
double calc_rho_n(double D, int idata, double *rho_n);
void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb);
double calc_rhoB(double xb, double yb, double zb);
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz);
int get_p_integral(int nji, double *ls, double *ks);
double getx2y(int n, double *x, double *y, double xin);
double getx2y_ist(int n, double *x, double *y, double xin, int *ist);
double getx2y_khi(int n, double *x, double *y, double xin, int *khi);
double interp_x(int n, double *F, double xst, double dx, double xreq);
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq);
void interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq);
}

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
