#include "genulens/model/coordinates.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"
namespace genulens {

void PopulationRuntime::initialize_mass_function(const InitialMassFunctionOptions &options) {
  nm = 1000;
  log_mass = (double*)calloc(nm + 1, sizeof(double *));
  mass_probability = (double*)calloc(nm + 1, sizeof(double *));
  mass_cumulative = (double*)calloc(nm + 1, sizeof(double *));
  mass_percentiles = (int*)calloc(22, sizeof(int *));
  store_IMF_nBs(1, log_mass, mass_probability, mass_cumulative, mass_percentiles,
                options.m0, options.m1, options.m2, options.m3, options.ml, options.mu,
                options.alpha1, options.alpha2, options.alpha3, options.alpha4, options.alpha0);
}

void PopulationRuntime::initialize_luminosity_functions(
    double source_i_min,
    double source_i_max,
    double source_vi_min,
    double source_vi_max,
    double ai_rc,
    double evi_rc) {
  const int narrow_lf_capacity = 960;
  if (source_i_max - source_i_min > 0 && source_vi_max - source_vi_min == 0 && ai_rc > 0) {
    CumuN_MIs = (double**)malloc(sizeof(double *) * ncomp);
    for (int i = 0; i < ncomp; i++) {
      CumuN_MIs[i] = (double*)calloc(narrow_lf_capacity, sizeof(double *));
    }
    MIs = (double*)malloc(sizeof(double *) * narrow_lf_capacity);
    nMIs = make_LFs(MIs, CumuN_MIs, log_mass, mass_cumulative);
    dILF = (MIs[nMIs - 1] - MIs[0]) / (nMIs - 1);
    has_luminosity_function = true;
  }

  if (source_i_max - source_i_min > 0 && source_vi_max - source_vi_min > 0 && ai_rc > 0 && evi_rc > 0) {
    double MIst = -5, MIen = 10, VIst = 0.0, VIen = 3.0;
    nMIs = 150;
    nVIs = 30;
    MIs = (double*)calloc(nMIs + 1, sizeof(double *));
    VIs = (double*)calloc(nVIs + 1, sizeof(double *));
    f_VI_Is = (double***)malloc(sizeof(double *) * ncomp);
    for (int i = 0; i < ncomp; i++) {
      f_VI_Is[i] = (double**)malloc(sizeof(double *) * (nVIs + 1));
      for (int j = 0; j < nVIs + 1; j++) {
        f_VI_Is[i][j] = (double*)calloc(nMIs + 1, sizeof(double *));
      }
    }
    store_VI_MI(MIst, MIen, nMIs, VIst, VIen, nVIs, MIs, VIs, f_VI_Is, log_mass, mass_cumulative);
    dILF = (double)(MIen - MIst) / nMIs;
    dVILF = (double)(VIen - VIst) / nVIs;
    has_color_magnitude_function = true;
  }
}

void PopulationRuntime::read_empirical_mass_luminosity() {
  int capacity = 60;
  empirical_masses = (double*)calloc(capacity, sizeof(double *));
  empirical_magnitudes = (double**)malloc(sizeof(double *) * empirical_columns);
  for (int i = 0; i < empirical_columns; i++) {
    empirical_magnitudes[i] = (double*)calloc(capacity, sizeof(double *));
  }
  char *file_MLemp = (char*)"input_files/MLemp.dat";
  empirical_count = read_MLemp(file_MLemp, empirical_masses, empirical_magnitudes);
}

void PopulationRuntime::release_luminosity_functions() {
  if (has_luminosity_function) {
    for (int i = 0; i < ncomp; i++) {
      free(CumuN_MIs[i]);
    }
    free(CumuN_MIs);
    free(MIs);
    has_luminosity_function = false;
  }
  if (has_color_magnitude_function) {
    for (int i = 0; i < ncomp; i++) {
      for (int j = 0; j < nVIs + 1; j++) {
        free(f_VI_Is[i][j]);
      }
      free(f_VI_Is[i]);
    }
    free(f_VI_Is);
    free(MIs);
    free(VIs);
    has_color_magnitude_function = false;
  }
}

void PopulationRuntime::release_all() {
  release_luminosity_functions();
  free(log_mass);
  free(mass_cumulative);
  free(mass_probability);
  free(mass_percentiles);
  if (empirical_magnitudes != nullptr) {
    for (int i = 0; i < empirical_columns; i++) {
      free(empirical_magnitudes[i]);
    }
    free(empirical_magnitudes);
  }
  free(empirical_masses);
}

void store_IMF_nBs(int B, double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles, double M0, double M1, double M2, double M3, double Ml, double Mu, double alpha1, double alpha2, double alpha3, double alpha4, double alpha0){
  /* Store IMF with a broken-power law form.
   * Update normalize factors for the density distribution if B == 1 */
  gmodel::IMFParameters imf_parameters;
  imf_parameters.m0 = M0;
  imf_parameters.m1 = M1;
  imf_parameters.m2 = M2;
  imf_parameters.m3 = M3;
  imf_parameters.ml = Ml;
  imf_parameters.mu = Mu;
  imf_parameters.alpha0 = alpha0;
  imf_parameters.alpha1 = alpha1;
  imf_parameters.alpha2 = alpha2;
  imf_parameters.alpha3 = alpha3;
  imf_parameters.alpha4 = alpha4;
  const auto mass_grid = gmodel::BrokenPowerLawIMF(imf_parameters).build_grid(nm, Ml, Mu);
  std::vector<double> PlogM_cum = mass_grid.cumulative_number;
  std::vector<double> PMlogM_cum = mass_grid.cumulative_mass;
  std::vector<double> PMlogM_cum_norm = mass_grid.cumulative_mass_norm;
  logMst = mass_grid.log_mass_start;
	dlogM = mass_grid.log_mass_step;
  for(int i=0;i<=nm;i++){
    logMass[i] = mass_grid.log_mass[i];
    PlogM[i] = mass_grid.probability_log_mass[i];
    PlogM_cum_norm[i] = mass_grid.cumulative_number_norm[i];
  }
  for(int i=0;i<22;i++){
    imptiles[i] = mass_grid.percentile_index[i];
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
    double PM   = interp_x(nm+1, PMlogM_cum_norm.data(), logMst, dlogM, logMdie);
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
  for (int i=0;i<9;i++){
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
    }else if (i == 7){  // Thick disk
      double logMdie = log10(MinidieD[nageD - 2]);
      double logMRG1 = log10(MRGstD[nageD - 2]);
      double logMRG2 = log10(MRGenD[nageD - 2]);
      double PM   = interp_x(nm+1, PMlogM_cum_norm.data(), logMst, dlogM, logMdie);
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
    }else{  // Stellar halo (added on 20241002)
      double logMdie = log10(MinidieD[nageD - 1]);
      double logMRG1 = log10(MRGstD[nageD - 1]);
      double logMRG2 = log10(MRGenD[nageD - 1]);
      double PM   = interp_x(nm+1, PMlogM_cum_norm.data(), logMst, dlogM, logMdie);
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
      double norm = rho0SHMS/PM;
      double rhoMS = rho0SHMS;
      double rhoWD = norm * PMWD;
      rho0SH = rhoMS + rhoWD;
      n0MSSH = rhoMS/aveMMS;
      double n0WD = rhoWD/aveMWD;
      n0SH   = n0MSSH + n0WD;
      n0RGSH = n0MSSH* PRG/P;
      // printf ("%d rho0= %.2e + %.2e = %.2e, n0= %.2e + %.2e = %.2e, n0RG= %.2e\n",i,rhoMS,rhoWD,rho0SH,n0MSSH,n0WD,n0SH,n0RGSH);
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
    double PM   = interp_x(nm+1, PMlogM_cum_norm.data(), logMst, dlogM, logMdie);
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
    double PM   = interp_x(nm+1, PMlogM_cum_norm.data(), logMst, dlogM, logMdie);
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
}

//----------------
void Mini2Mrem (double *pout, double Mini, int mean) {  // mean = 1: give mean, 0: give random
  const auto remnant = gmodel::RemnantMassModel(MiniWDmax).evolve(Mini, mean == 1, *active_state->runtime.rng);
  pout[0] = remnant.mass_msun;
  pout[1] = remnant.remnant_type;
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
int make_LFs(double *MIs_arg, double **CumuN_MIs_arg, double *logMass, double *PlogM_cum_norm)
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
             (icomp == 9) ? "input_files/iso-NSD.dat" :
                            "input_files/iso-halo.dat";
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
            (icomp == 9) ? 6  :
                           2;
     narry = (icomp == 0) ? 465 :
             (icomp == 1) ? 581 :
             (icomp == 2) ? 1885 :
             (icomp == 3) ? 552 :
             (icomp == 4) ? 362 :
             (icomp == 5) ? 323 :
             (icomp == 6) ? 296 :
             (icomp == 7) ? 197 :
             (icomp == 8) ? 766 :
             (icomp == 9) ? 340 :
                            190;
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
            (icomp == 9) ? 100 :
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
                    : 2; // Thick or halo
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
void store_VI_MI(double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI, double *MIs_arg, double *VIs_arg, double ***f_VI_Is_arg, double *logMass, double *PlogM_cum_norm){
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
            (icomp == 9) ? "input_files/iso-NSD.dat" :
                           "input_files/iso-halo.dat";
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
           (icomp == 9) ? 6  :
                          2;
    narry = (icomp == 0) ? 465 :
            (icomp == 1) ? 581 :
            (icomp == 2) ? 1885 :
            (icomp == 3) ? 552 :
            (icomp == 4) ? 362 :
            (icomp == 5) ? 323 :
            (icomp == 6) ? 296 :
            (icomp == 7) ? 197 :
            (icomp == 8) ? 766 :
            (icomp == 9) ? 340 :
                           190;
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
                   : 2; // Thick or Halo
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

} // namespace genulens
