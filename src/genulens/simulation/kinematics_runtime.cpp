#include "genulens/model/coordinates.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"
namespace genulens {

void KinematicRuntimeTables::initialize_shu_distribution() {
  n_z = (zenShu - zstShu) / dzShu + 1;
  n_R = (RenShu - RstShu) / dRShu + 1;
  fgsShu = (double****)malloc(sizeof(double *) * n_z);
  PRRgShus = (double****)malloc(sizeof(double *) * n_z);
  cumu_PRRgs = (double****)malloc(sizeof(double *) * n_z);
  n_fgsShu = (int***)malloc(sizeof(int *) * n_z);
  kptiles = (int****)malloc(sizeof(int *) * n_z);
  for (int i = 0; i < n_z; i++) {
    fgsShu[i] = (double***)malloc(sizeof(double *) * n_R);
    PRRgShus[i] = (double***)malloc(sizeof(double *) * n_R);
    cumu_PRRgs[i] = (double***)malloc(sizeof(double *) * n_R);
    n_fgsShu[i] = (int**)malloc(sizeof(int *) * n_R);
    kptiles[i] = (int***)malloc(sizeof(int *) * n_R);
    for (int j = 0; j < n_R; j++) {
      fgsShu[i][j] = (double**)malloc(sizeof(double *) * n_disk);
      PRRgShus[i][j] = (double**)malloc(sizeof(double *) * n_disk);
      cumu_PRRgs[i][j] = (double**)malloc(sizeof(double *) * n_disk);
      kptiles[i][j] = (int**)malloc(sizeof(int *) * n_disk);
      n_fgsShu[i][j] = (int*)calloc(n_disk, sizeof(int *));
      for (int k = 0; k < n_disk; k++) {
        fgsShu[i][j][k] = (double*)calloc(n_fg, sizeof(double *));
        PRRgShus[i][j][k] = (double*)calloc(n_fg, sizeof(double *));
        cumu_PRRgs[i][j][k] = (double*)calloc(n_fg, sizeof(double *));
        kptiles[i][j][k] = (int*)calloc(22, sizeof(int *));
      }
    }
  }
  char *fileVc = (char*)"input_files/Rotcurve_BG16.dat";
  store_cumuP_Shu(fileVc);
}

void KinematicRuntimeTables::release_all() {
  for (int i = 0; i < n_z; i++) {
    for (int j = 0; j < n_R; j++) {
      for (int k = 0; k < n_disk; k++) {
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
}

void NsdMomentRuntime::initialize_if_enabled() {
  nzND = (zenND - zstND) / dzND + 1.5;
  nRND = (RenND - RstND) / dRND + 1.5;
  if (ND != 3) {
    return;
  }
  logrhoNDs = (double**)malloc(sizeof(double *) * nzND);
  vphiNDs = (double**)malloc(sizeof(double *) * nzND);
  corRzNDs = (double**)malloc(sizeof(double *) * nzND);
  logsigvNDs = (double***)malloc(sizeof(double *) * nzND);
  for (int i = 0; i < nzND; i++) {
    logrhoNDs[i] = (double*)calloc(nRND, sizeof(double *));
    vphiNDs[i] = (double*)calloc(nRND, sizeof(double *));
    corRzNDs[i] = (double*)calloc(nRND, sizeof(double *));
    logsigvNDs[i] = (double**)malloc(sizeof(double *) * nRND);
    for (int j = 0; j < nRND; j++) {
      logsigvNDs[i][j] = (double*)calloc(3, sizeof(double *));
    }
  }
  char *fileND = (char*)"input_files/NSD_moments.dat";
  store_NSDmoments(fileND);
}

void NsdMomentRuntime::release_if_enabled() {
  if (ND != 3) {
    return;
  }
  for (int i = 0; i < nzND; i++) {
    for (int j = 0; j < nRND; j++) {
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
    double facVcz = 1 + 0.0374*pow(0.001*abs(z), 1.34); // Eq. (22) of Sharma et al. 2014, ApJ, 793, 51
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
  if (i < 8 && B14disk == 1){
    double aveV = (i < 7) ? 218.0 : 170.0;
    double sigU = (i < 7) ?  39.9 :  67.0;
    double sigV = (i < 7) ?  27.9 :  51.0;
    double sigW = (i < 7) ?  19.1 :  42.0;
    do{
      double vR   =    0 + gasdev()*sigU; // radial velocity
      double vphi = aveV + gasdev()*sigV;
      vx = -vphi * y/R + vR * x/R; // x/R = cosphi, y/R = sinphi
      vy =  vphi * x/R + vR * y/R; // x/R = cosphi, y/R = sinphi
      vz =    0 + gasdev()*sigW; // vertical velocity
    }while (vx*vx + vy*vy + vz*vz > vescd*vescd);
  }else if (i < 8){
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
  }else if (i == 10){ // Stellar halo
    double vR   = sigU_SH * gasdev();
    double vphi = sigV_SH * gasdev();
    vx = -vphi * y/R + vR * x/R; // x/R = cosphi, y/R = sinphi
    vy =  vphi * x/R + vR * y/R; // x/R = cosphi, y/R = sinphi
    vz =    0 + gasdev()*sigW_SH; // vertical velocity
  }else{ // bar & NSD (when ND <= 2)
    double vrot = 0.001 * Omega_p * R; // km/s/kpc -> km/s/pc
    if (B14vbar == 1 && vrot > 90)  vrot = 90;
    double xb =  x * costheta + y * sintheta;
    double yb = -x * sintheta + y * costheta;
    double zb =  z;                          
    double sigvbs[3] = {}, sigx, sigy, sigz;
    if (B14vbar==1){
      double Mvb_xrot =  -vrot * y/R;
      double Mvb_yrot =   vrot * x/R;
      sigx = sqrt(sigx_vb*sigx_vb - Mvb_xrot*Mvb_xrot);
      sigy = sqrt(sigy_vb*sigy_vb - Mvb_yrot*Mvb_yrot);
      sigz = sigz_vb;
    }else{
      void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
      calc_sigvb(xb, yb, zb, sigvbs);
      sigx = sqrt(sigvbs[0]*sigvbs[0] * costheta*costheta + sigvbs[1]*sigvbs[1] * sintheta*sintheta);
      sigy = sqrt(sigvbs[0]*sigvbs[0] * sintheta*sintheta + sigvbs[1]*sigvbs[1] * costheta*costheta);
      sigz = sigvbs[2];
    }
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
   const auto separation = gmodel::BinaryLensSampler().sample_projected_separation(M1, M2, coeff, *active_state->runtime.rng);
   pout[0] = separation.log_separation_au;
   pout[1] = separation.projected_separation_au;
}

double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv){ 
  return genulens::math::Interpolation::inverse_cumulative_linear_density(n, x, F, f, Freq, ist, inv);
}
double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd){ 
  double getx2y(int n, double *x, double *y, double xin);
  if (fg <= 0) return 0;
  double  calc_faca(double Rg, double hsigU, int rd, double a0);
  double calc_SigRg(double Rg, double hsigU, int rd, double a0);
  double    calc_gc(double c);
  double Rg = R*fg;
  double vc = getx2y(nVcs, Rcs, Vcs, Rg) / (1 + 0.0374*pow(0.001*abs(z), 1.34));
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

} // namespace genulens
