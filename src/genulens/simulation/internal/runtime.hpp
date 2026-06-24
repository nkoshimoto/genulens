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

  void initialize_mass_function(RunContext &ctx, const InitialMassFunctionOptions &options);
  void initialize_luminosity_functions(
      RunContext &ctx,
      double source_i_min,
      double source_i_max,
      double source_vi_min,
      double source_vi_max,
      double ai_rc,
      double evi_rc);
  void read_empirical_mass_luminosity();
  void release_luminosity_functions(RunContext &ctx);
  void release_all(RunContext &ctx);
};

struct KinematicRuntimeTables {
  int n_z = 0;
  int n_R = 0;
  int n_disk = 8;
  int n_fg = 100;

  void initialize_shu_distribution(RunContext &ctx);
  void release_all(RunContext &ctx);
};

struct NsdMomentRuntime {
  void initialize_if_enabled(RunContext &ctx);
  void release_if_enabled(RunContext &ctx);
};

void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA);
void calc_opticaldepth(RunContext &ctx, double *tauall, double *Nsall, int idata, int Dsmax21, double AI0, double hscale, double Isst, double Isen);
void store_NSDmoments(RunContext &ctx, char *infile);
double like_obs(RunContext &ctx, double mod, double obs, double err, double fe, int det, int UNIFORM);
void store_IMF_nBs(RunContext &ctx, int B, double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles, double M0, double M1, double M2, double M3, double Mbr, double Ml, double Mu, double alpha1, double alpha2, double alpha3, double alpha4, double alpha5, double alpha0);
void Mini2Mrem(RunContext &ctx, double *pout, double Mini, int mean);
double fLF_detect(RunContext &ctx, double extI, double Imin, double Imax, int idisk);
double fIVI_detect(RunContext &ctx, double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk);
void store_cumuP_Shu(RunContext &ctx, char *infile);
void get_PRRGmax2(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
void calc_dpdfg(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd);
void get_vxyz_ran(RunContext &ctx, double *vxyz, int i, double tau, double D, double lD, double bD);
void getaproj(RunContext &ctx, double *pout, double M1, double M2, int coeff);
double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
int read_MLemp(char *infile, double *M_emps, double **Mag_emps);
int make_LFs(RunContext &ctx, double *MIs_arg, double **CumuN_MIs_arg, double *logMass, double *PlogM_cum_norm);
void store_VI_MI(RunContext &ctx, double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI, double *MIs_arg, double *VIs_arg, double ***f_VI_Is_arg, double *logMass, double *PlogM_cum_norm);
double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd);
double calc_gc(double c);
double calc_SigRg(double Rg, double hsigU, int rd, double a0);
double calc_faca(double Rg, double hsigU, int rd, double a0);
double crude_integrate(RunContext &ctx, double xmax, double ymax, double zmax, int nbun);
void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
void cross(double *c, double *a, double *b);
double dot(double *a, double *b);
void norm_vec(double *a);
double calc_rho_n(RunContext &ctx, double D, int idata, double *rho_n);
void calc_rho_each(RunContext &ctx, double D, int idata, double *rhos, double *xyz, double *xyb);
double calc_rhoB(double xb, double yb, double zb);
void Dlb2xyz(const RunContext &ctx, double D, double lD, double bD, double Rsun, double *xyz);
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
#define NDATAMAX 8000000000LL // take ~6hours?
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
