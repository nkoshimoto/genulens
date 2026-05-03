#include "genulens/simulation/internal/runtime.hpp"

namespace genulens {

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
  // rho_n[1] += n0MSb * rhos[8] + n0ND   * rhos[9];
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

  // Stellar Halo (Spheroid from Robin+03) added on 20241002
  double zSH  = z/epsSH;
  double aSH2 = R*R + zSH*zSH;
  if (aSH2 < acSH2) aSH2 = acSH2; // a = ac = 500 pc
  double aSHsun_inv = R0/sqrt(aSH2); // since rho0SHMS is the value at Solar position
  rhos[10] = pow(aSHsun_inv, alphaSH); // a^-alpha = (1/a)^alpha
  
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

} // namespace genulens
