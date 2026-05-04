#include "genulens/model/coordinates.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"
namespace genulens {

void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA){
  const auto pa = gmodel::CoordinateTransformer::position_angle(gl, gb);
  *PA = pa.degrees;
  *cosPA = pa.cos_pa;
  *sinPA = pa.sin_pa;
}

//----------------
void calc_opticaldepth(RunContext &ctx, double *tauall, double *Nsall, int idata, int Dsmax21, double AI0, double hscale, double Isst, double Isen)
{
  active_state = &ctx;
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
  double calc_rho_n(RunContext &ctx, double D, int idata, double *rho_n); // return rho(D) & n(D)
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
          wtDB = calc_rho_n(ctx, Dl, idata, rho_n);
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
          wtDB = calc_rho_n(ctx, Dl, idata, rho_n);
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
          wtDB = calc_rho_n(ctx, Dl, idata, rho_n);
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
  double fLF_detect(RunContext &ctx, double extI, double Imin, double Imax, int idisk);
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
     double fMagD =  (1 - wtD) * fLF_detect(ctx, extI, Isst, Isen, (int)floor(m_idisk))
                   +      wtD  * fLF_detect(ctx, extI, Isst, Isen, (int) ceil(m_idisk));

     // Bulge/bar fraction of 14 < I < 21
     // Assume same fMag for NSD as Bulge
     double fMagB =  fLF_detect(ctx, extI, Isst, Isen, 8);

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
     double fMagD =  (1 - wtD) * fLF_detect(ctx, extI, Isst, Isen, (int)floor(m_idisk))
                   +      wtD  * fLF_detect(ctx, extI, Isst, Isen, (int) ceil(m_idisk));

     // Bulge/bar fraction of 14 < I < 21
     // Assume same fMag for NSD as Bulge
     double fMagB =  fLF_detect(ctx, extI, Isst, Isen, 8);

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
void store_NSDmoments(RunContext &ctx, char *infile) // Read input_files/NSD_moments.dat
{
  active_state = &ctx;
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
double like_obs(RunContext &ctx, double mod, double obs, double err, double fe, int det, int UNIFORM){
  genulens::ObservationConstraint constraint;
  constraint.observed = obs;
  constraint.error = err;
  constraint.error_inflation = fe;
  constraint.detection_mode = static_cast<genulens::DetectionMode>(det);
  constraint.uniform = (UNIFORM == 1);
  return genulens::ObservationLikelihood(constraint).accept_weight(mod, *ctx.runtime.rng);
}

} // namespace genulens
