/* galactic_model.cpp
 * Shared Galactic model functions for pre_gapmoe tools.
 * Extracted from genulens.cpp (Koshimoto, Baba & Bennett 2021).
 * Replaces Numerical Recipes ran1/gasdev with GSL equivalents.
 */
#include "galactic_model.h"
#include "option.h"

// =========================================================
// Global variable definitions
// =========================================================

const gsl_rng_type *T_rng;
gsl_rng *r_rng;

int ncomp = 11;

int DISK, hDISK, addX, model;
int ND = 0, SH = 1;

double R0, thetaD;
double x0_1, y0_1, z0_1 = 0;
double C1, C2, C3;
double Rc, frho0b, zb_c;
double costheta, sintheta;

double x0_X, y0_X, z0_X = 0;
double C1_X, C2_X;
double b_zX, b_zY = 0, fX, Rc_X;

int    Rd[3]   = {5000, 2600, 2200};
int    Rh = 3740, Rdbreak = 5300, nh = 1;
double zd[8]   = {61.47, 141.84, 224.26, 292.36, 372.85, 440.71, 445.37, 903.12};
double zd45[8] = {36.88,  85.10, 134.55, 175.41, 223.71, 264.42, 267.22, 903.12};
double y0d[3];

double tSFR = 7.0;
double medtauds[8] = {0.075273, 0.586449, 1.516357, 2.516884, 4.068387, 6.069263, 8.656024, 12};

double hsigUt, hsigWt, hsigUT, hsigWT;
double betaU, betaW, sigU10d, sigW10d, sigU0td, sigW0td;

int    model_vb, model_vbz;
double Omega_p;
double x0_vb, y0_vb, z0_vb, C1_vb, C2_vb, C3_vb;
double sigx_vb, sigy_vb, sigz_vb;
double sigx_vb0, sigy_vb0, sigz_vb0;
double vx_str, y0_str;
double x0_vbz, y0_vbz, z0_vbz, C1_vbz, C2_vbz, C3_vbz;

double vxsun = -10.0, Vsun = 11.0, vzsun = 7.0, vysun = 243.0;

double rhot0;
double rho0d[8] = {5.16e-03+3.10e-04, 5.00e-03+5.09e-04, 3.85e-03+5.42e-04, 3.18e-03+5.54e-04,
                   5.84e-03+1.21e-03,  6.24e-03+1.51e-03,  1.27e-02+3.49e-03,  1.68e-03+6.02e-04};
double n0d[8]   = {1.51e-02+1.12e-04, 1.66e-02+3.22e-04, 1.40e-02+4.39e-04, 1.22e-02+5.15e-04,
                   2.36e-02+1.25e-03,  2.63e-02+1.67e-03,  5.55e-02+4.08e-03,  7.91e-03+7.81e-04};
double n0MSd[8] = {1.51e-02, 1.66e-02, 1.40e-02, 1.22e-02, 2.36e-02, 2.63e-02, 5.55e-02, 7.91e-03};
double n0RGd[8] = {7.09e-06, 3.40e-05, 4.32e-05, 2.16e-05, 6.60e-05, 6.19e-05, 1.29e-04, 9.38e-06};

double rho0b, n0MSb, n0RGb, n0b;
double fb_MS  = 1.62/2.07;
double m2nb_MS = 1/0.227943;
double m2nb_WD = 1/0.847318;
double nMS2nRGb = 2.33232e-03;

int    x0ND = 250, y0ND = 125, z0ND = 50;
double C1ND = 2, rho0ND, n0MSND, n0RGND, n0ND;
double fND_MS = 0, m2nND_MS = 0, m2nND_WD = 0, nMS2nRGND = 0;
double **logrhoNDs, **vphiNDs, ***logsigvNDs, **corRzNDs;
double zstND = 0, zenND = 400, dzND = 5;
double RstND = 0, RenND = 1000, dRND = 5;
int    nzND, nRND;

double rho0SHMS = 9.32e-06;
double epsSH    = 0.76;
double alphaSH  = 2.44;
double acSH2    = 500*500;
double rho0SH, n0MSSH, n0RGSH, n0SH;
double sigU_SH, sigV_SH, sigW_SH;

int    agesD[250], agesB[50], agesND[10];
double MinidieD[250], MinidieB[50], MinidieND[10];
int    nageD = 0, nageB = 0, nageND = 0;
double MiniWDmax = 9;
double mageB = 9, sageB = 1, mageND = 7, sageND = 1;

int    nm;
double logMst, dlogM;

double ****fgsShu    = NULL;
double ****PRRgShus  = NULL;
double ****cumu_PRRgs= NULL;
int   ***n_fgsShu    = NULL;
int   ****kptiles    = NULL;
int    zstShu =   0, zenShu = 3600, dzShu = 200;
int    RstShu = 500, RenShu = 12200, dRShu = 100;

int    nVcs = 0;
double Rcs[60], Vcs[60];

double *lDs, *bDs;
double xyzSgrA[3] = {};

int B14disk = 0, B14vbar = 0;

// IMF table kept alive after init_galactic_model for use by calc_mass_dist
double *g_logMass        = NULL;
double *g_PlogM          = NULL;
double *g_PlogM_cum_norm = NULL;

// =========================================================
// GSL random number wrappers
// =========================================================

double ran1() {
    return gsl_rng_uniform(r_rng);
}

double gasdev() {
    return gsl_ran_ugaussian(r_rng);
}

// =========================================================
// Utility: linear interpolation
// =========================================================

double getx2y(int n, double *x, double *y, double xin)
{
    double xmin, xmax;
    if (x[0] < x[n-1]) { xmin = x[0]; xmax = x[n-1]; }
    else                { xmin = x[n-1]; xmax = x[0]; }
    if (xmin > xin) return 0;
    if (xmax < xin) return 0;
    for (int i = 1; i < n; i++) {
        if ((x[i] <= xin && x[i-1] >= xin) || (x[i] >= xin && x[i-1] <= xin)) {
            return (y[i]-y[i-1])/(x[i]-x[i-1])*(xin-x[i-1]) + y[i-1];
        }
    }
    return 0;
}

double getx2y_ist(int n, double *x, double *y, double xin, int *ist)
{
    for (int i = *ist; i < n; i++) {
        if (i == 0) continue;
        if ((x[i] <= xin && x[i-1] >= xin) || (x[i] >= xin && x[i-1] <= xin)) {
            double yreq = (y[i]-y[i-1])/(x[i]-x[i-1])*(xin-x[i-1]) + y[i-1];
            *ist = i;
            return yreq;
        }
    }
    return 0;
}

double getx2y_khi(int n, double *x, double *y, double xin, int *khi)
{
    double xmin, xmax;
    if (x[0] < x[n-1]) { xmin = x[0]; xmax = x[n-1]; }
    else                { xmin = x[n-1]; xmax = x[0]; }
    if (xmin > xin) return 0;
    if (xmax < xin) return 0;
    int klo;
    if (*khi > 0) {
        klo = *khi - 1;
    } else {
        klo = 0;
        *khi = n-1;
    }
    while (*khi - klo > 1) {
        int k = (*khi + klo) >> 1;
        if (x[k] > xin) *khi = k;
        else              klo = k;
    }
    if (x[*khi] == x[klo]) return y[klo];
    return (y[*khi]-y[klo])/(x[*khi]-x[klo])*(xin-x[klo]) + y[klo];
}

double interp_x(int n, double *F, double xst, double dx, double xreq)
{
    int    ix   = (int)((xreq - xst)/dx);
    double xres = (xreq - xst)/dx - ix;
    if (ix < 0 || ix > n-1) return 0;
    if (ix+1 > n-1) return F[ix];
    return F[ix]*(1-xres) + F[ix+1]*xres;
}

double interp_xy(int nx, int ny, double **F,
                 double xst, double yst, double dx, double dy,
                 double xreq, double yreq)
{
    int    ix   = (int)((xreq - xst)/dx);
    double xres = (xreq - xst)/dx - ix;
    int    iy   = (int)((yreq - yst)/dy);
    double yres = (yreq - yst)/dy - iy;
    if (ix < 0 || ix > nx-1 || iy < 0 || iy > ny-1) return 0;
    if (ix+1 > nx-1 && iy+1 > ny-1) return F[ix][iy];
    if (ix+1 > nx-1) return F[ix][iy]*(1-yres) + F[ix][iy+1]*yres;
    if (iy+1 > ny-1) return F[ix][iy]*(1-xres) + F[ix+1][iy]*xres;
    return (1-xres)*(1-yres)*F[ix][iy]   + xres*(1-yres)*F[ix+1][iy]
         + (1-xres)*    yres*F[ix][iy+1] + xres*    yres*F[ix+1][iy+1];
}

void interp_xy_coeff(int nx, int ny, double *as,
                     double xst, double yst, double dx, double dy,
                     double xreq, double yreq)
{
    int    ix   = (int)((xreq - xst)/dx);
    double xres = (xreq - xst)/dx - ix;
    int    iy   = (int)((yreq - yst)/dy);
    double yres = (yreq - yst)/dy - iy;
    if (ix < 0 || ix > nx-1 || iy < 0 || iy > ny-1) {
        as[0] = as[1] = as[2] = as[3] = 0;
    } else if (ix+1 > nx-1 && iy+1 > ny-1) {
        as[0] = 1; as[1] = as[2] = as[3] = 0;
    } else if (ix+1 > nx-1) {
        as[0] = 1-yres; as[2] = yres; as[1] = as[3] = 0;
    } else if (iy+1 > ny-1) {
        as[0] = 1-xres; as[1] = xres; as[2] = as[3] = 0;
    } else {
        as[0] = (1-xres)*(1-yres);
        as[1] =    xres *(1-yres);
        as[2] = (1-xres)*   yres ;
        as[3] =    xres *   yres ;
    }
}

// =========================================================
// Linear algebra helpers
// =========================================================

void cross(double *c, double *a, double *b) {
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
}

double dot(double *a, double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void norm_vec(double *a) {
    double norm = dot(a, a);
    for (int i = 0; i < 3; i++) a[i] /= sqrt(norm);
}

// =========================================================
// Coordinate conversion
// =========================================================

void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz)
{
    double cosbsun = cos(zsun/Rsun), sinbsun = sin(zsun/Rsun);
    double cosb = cos(bD/180.0*PI), sinb = sin(bD/180.0*PI);
    double cosl = cos(lD/180.0*PI), sinl = sin(lD/180.0*PI);
    double xtmp = Rsun - D*cosb*cosl;
    double ytmp = D*cosb*sinl;
    double ztmp = D*sinb;
    xyz[0] = xtmp - xyzSgrA[0];
    xyz[1] = ytmp - xyzSgrA[1];
    xyz[2] = ztmp*cosbsun + xtmp*sinbsun - xyzSgrA[2];
}

void calc_PA(double gl, double gb, double *PA, double *cosPA, double *sinPA)
{
    double lNP = 122.9320, bNP = 27.1284;
    double cosl   = cos(gl/180.*PI), sinl   = sin(gl/180.*PI);
    double cosb   = cos(gb/180.*PI), sinb   = sin(gb/180.*PI);
    double coslNP = cos(lNP/180.*PI), sinlNP = sin(lNP/180.*PI);
    double cosbNP = cos(bNP/180.*PI), sinbNP = sin(bNP/180.*PI);
    double nvector[3] = {cosb*cosl, cosb*sinl, sinb};
    double galNP[3]   = {0, 0, 1};
    double eqNP[3]    = {cosbNP*coslNP, cosbNP*sinlNP, sinbNP};
    double elvector[3] = {}, ebvector[3] = {}, eEvector[3] = {}, eNvector[3] = {};
    cross(elvector, galNP, nvector);
    norm_vec(elvector);
    cross(ebvector, nvector, elvector);
    cross(eEvector, eqNP, nvector);
    norm_vec(eEvector);
    cross(eNvector, nvector, eEvector);
    double crosstmp[3] = {};
    *cosPA = dot(elvector, eEvector);
    cross(crosstmp, elvector, eEvector);
    *sinPA = -dot(nvector, crosstmp);
    *PA    = 180*atan2(*sinPA, *cosPA)/PI;
}

// =========================================================
// Density: bulge
// =========================================================

double calc_rhoB(double xb, double yb, double zb)
{
    double xn, yn, zn, R, Rs, rs, rho = 0, rhoX = 0;
    R = sqrt(xb*xb + yb*yb);
    if (model >= 4 && model <= 8) {
        xn = fabs(xb/x0_1), yn = fabs(yb/y0_1), zn = fabs(zb/z0_1);
        Rs = pow(pow(xn,C1) + pow(yn,C1), 1.0/C1);
        rs = pow(pow(Rs,C2) + pow(zn,C2), 1.0/C2);
        if (rs == 0 && model == 8) rs = 0.0001;
        rho = (model == 5) ? exp(-rs)
            : (model == 6) ? exp(-0.5*rs*rs)
            : (model == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2)
            : (model == 4) ? exp(-pow(rs, C3))
            : 0;
    }
    if (R  >= Rc)   rho *= exp(-0.5*(R-Rc)*(R-Rc)/srob/srob);
    if (fabs(zb) >= zb_c) rho *= exp(-0.5*(fabs(zb)-zb_c)*(fabs(zb)-zb_c)/200.0/200.0);

    if (addX >= 5) {
        xn = fabs((xb-b_zX*zb)/x0_X), yn = fabs((yb-b_zY*zb)/y0_X), zn = fabs(zb/z0_X);
        rs = pow(pow(pow(xn,C1_X)+pow(yn,C1_X), C2_X/C1_X) + pow(zn,C2_X), 1.0/C2_X);
        rhoX  = (addX == 5) ? exp(-rs)
              : (addX == 6) ? exp(-0.5*rs*rs)
              : (addX == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2) : 0;
        xn = fabs((xb+b_zX*zb)/x0_X), yn = fabs((yb-b_zY*zb)/y0_X);
        rs = pow(pow(pow(xn,C1_X)+pow(yn,C1_X), C2_X/C1_X) + pow(zn,C2_X), 1.0/C2_X);
        rhoX += (addX == 5) ? exp(-rs)
              : (addX == 6) ? exp(-0.5*rs*rs)
              : (addX == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2) : 0;
        if (b_zY > 0.0) {
            xn = fabs((xb-b_zX*zb)/x0_X), yn = fabs((yb+b_zY*zb)/y0_X);
            rs = pow(pow(pow(xn,C1_X)+pow(yn,C1_X), C2_X/C1_X) + pow(zn,C2_X), 1.0/C2_X);
            rhoX += (addX == 5) ? exp(-rs)
                  : (addX == 6) ? exp(-0.5*rs*rs)
                  : (addX == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2) : 0;
            xn = fabs((xb+b_zX*zb)/x0_X), yn = fabs((yb+b_zY*zb)/y0_X);
            rs = pow(pow(pow(xn,C1_X)+pow(yn,C1_X), C2_X/C1_X) + pow(zn,C2_X), 1.0/C2_X);
            rhoX += (addX == 5) ? exp(-rs)
                  : (addX == 6) ? exp(-0.5*rs*rs)
                  : (addX == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2) : 0;
        }
        rhoX *= fX;
        if (R >= Rc_X) rhoX *= exp(-0.5*(R-Rc_X)*(R-Rc_X)/srob/srob);
        rho += rhoX;
    }
    return rho;
}

// =========================================================
// Density: all components
// Returns shape function rhos[i], normalized so rhos[i]=1 at solar position
// Multiply by n0MSd[i] / n0d[i] / etc. to get physical density
// =========================================================

void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb)
{
    double lD = lDs[idata], bD = bDs[idata];
    Dlb2xyz(D, lD, bD, R0, xyz);
    double x = xyz[0], y = xyz[1], z = xyz[2];
    double R = sqrt(x*x + y*y);
    for (int i = 0; i < ncomp; i++) rhos[i] = 0;

    // Disk (components 0-7)
    if (DISK > 0) {
        for (int idisk = 0; idisk < 8; idisk++) {
            double zdtmp = (hDISK == 0) ? zd[idisk]
                : (R > 4500) ? zd[idisk] + (R-R0)*(zd[idisk]-zd45[idisk])/(R0-4500)
                             : zd45[idisk];
            double rhotmp = (idisk < 7) ? 4.0/(exp(2*z/zdtmp)+exp(-2*z/zdtmp)+2)
                                        : exp(-fabs(z)/zd[idisk]);
            int itmp = (idisk == 0) ? 0 : (idisk < 7) ? 1 : 2;
            rhotmp *= zd[idisk]/zdtmp;
            if (DISK == 1) rhotmp *= exp(-R/Rd[itmp] - pow((double)Rh/R, nh));
            if (DISK == 2) rhotmp = (R > Rdbreak) ? rhotmp * exp(-R/Rd[itmp])
                                                   : rhotmp * exp(-(double)Rdbreak/Rd[itmp]);
            if (DISK == 3) rhotmp *= exp(-R/Rd[itmp]);
            rhos[idisk] = rhotmp / y0d[itmp];
        }
    }

    // Bar (component 8)
    double xb =  x*costheta + y*sintheta;
    double yb = -x*sintheta + y*costheta;
    double zb =  z;
    rhos[8] = calc_rhoB(xb, yb, zb);
    xyb[0] = xb; xyb[1] = yb;

    // NSD (component 9)
    if (ND > 0) {
        if (ND == 3) {
            if (R <= RenND-30 && fabs(z) <= zenND-20) {
                rhos[9] = pow(10.0, interp_xy(nzND, nRND, logrhoNDs, zstND, RstND, dzND, dRND, fabs(z), R));
            } else {
                rhos[9] = 0;
            }
        } else {
            double xn = fabs(xb/x0ND), yn = fabs(yb/y0ND), zn = fabs(zb/z0ND);
            double rs = pow(pow(xn,C1ND)+pow(yn,C1ND), 1.0/C1ND) + zn;
            rhos[9] = exp(-rs);
        }
    }

    // Stellar halo (component 10, Robin+03)
    double zSH  = z/epsSH;
    double aSH2 = R*R + zSH*zSH;
    if (aSH2 < acSH2) aSH2 = acSH2;
    rhos[10] = pow(R0/sqrt(aSH2), alphaSH);
}

// =========================================================
// Integration
// =========================================================

int get_p_integral(int nji, double *ls, double *ks)
{
    int nmin, i = 0, j = 0;
    nji = (nji <= 1) ? 1 : (nji <= 2) ? 2 : (nji <= 4) ? 4
        : (nji <= 6) ? 6 : (nji <= 8) ? 8 : 10;
    if (nji == 1) {
        ls[i++] = 0.0; ks[j++] = 0.5; nmin = 1;
    }
    if (nji == 2) {
        ls[i++]=0.0; ls[i++]=0.5; ls[i++]=1.0;
        ks[j++]=3.0/12; ks[j++]=4.0/12; ks[j++]=11.0/12;
        nmin = 3;
    }
    if (nji == 4) {
        ls[i++]=0.0;
        ls[i++]=1./4; ls[i++]=1./2; ls[i++]=3./4; ls[i++]=1.0;
        ls[i++]=3./2; ls[i++]=2.0; ls[i++]=9./4; ls[i++]=3.0;
        ks[j++]=70./360;
        ks[j++]=32./360; ks[j++]=76./360; ks[j++]=128./360; ks[j++]=187./360;
        ks[j++]=100./360; ks[j++]=218./360; ks[j++]=96./360; ks[j++]=353./360;
        nmin = 7;
    }
    // Higher-order formulas omitted for brevity; nji=2 (Simpson) is sufficient for normalization
    return nmin;
}

double crude_integrate(double xmax, double ymax, double zmax, int nbun)
{
    double *ls, *ks;
    int nji = 2;
    int narry = 3;
    ls = (double*)malloc(sizeof(double)*narry);
    ks = (double*)malloc(sizeof(double)*narry);
    int nmin = get_p_integral(nji, ls, ks);
    if (nbun < nmin) nbun = nmin;
    int ncalc = nbun + 1 + 2*narry - 2*nji;
    double *rhosumz  = (double*)malloc(sizeof(double)*ncalc);
    double *rhosumyz = (double*)malloc(sizeof(double)*ncalc);
    double dx = xmax/nbun, dy = ymax/nbun, dz = zmax/nbun;
    double totalmass = 0;
    for (int ix = 0; ix < ncalc; ix++) {
        int ixtmp = ix - 2*narry + nji;
        double xb = (ix >= 2*narry) ? dx*ixtmp
                  : (ix%2==0) ? dx*ls[ix/2] : xmax - dx*ls[ix/2];
        for (int iy = 0; iy < ncalc; iy++) {
            int iytmp = iy - 2*narry + nji;
            double yb = (iy >= 2*narry) ? dy*iytmp
                      : (iy%2==0) ? dy*ls[iy/2] : ymax - dy*ls[iy/2];
            rhosumz[iy] = 0;
            for (int j = 0; j < narry; j++) {
                double dztmp = dz*ls[j];
                rhosumz[iy] += (calc_rhoB(xb,yb,dztmp) + calc_rhoB(xb,yb,zmax-dztmp)) * ks[j];
            }
            for (int j = nji; j <= nbun-nji; j++) {
                rhosumz[iy] += calc_rhoB(xb, yb, dz*j);
            }
            rhosumz[iy] *= dz;
        }
        rhosumyz[ix] = 0;
        for (int j = 0; j < narry; j++) {
            rhosumyz[ix] += (rhosumz[2*j] + rhosumz[2*j+1]) * ks[j];
        }
        for (int j = 2*narry; j < ncalc; j++) rhosumyz[ix] += rhosumz[j];
        rhosumyz[ix] *= dy;
    }
    for (int j = 0; j < narry; j++) {
        totalmass += (rhosumyz[2*j] + rhosumyz[2*j+1]) * ks[j];
    }
    for (int j = 2*narry; j < ncalc; j++) totalmass += rhosumyz[j];
    totalmass *= dx * 8; // factor 8 for 8 octants
    free(ls); free(ks); free(rhosumz); free(rhosumyz);
    return totalmass;
}

// =========================================================
// NSD moment data loader
// =========================================================

void store_NSDmoments(char *infile)
{
    FILE *fp;
    char line[1000];
    char *words[100];
    if ((fp = fopen(infile, "r")) == NULL) {
        printf("can't open %s\n", infile);
        exit(1);
    }
    int iRz = 0;
    while (fgets(line, 1000, fp) != NULL) {
        split((char*)" ", line, words);
        if (*words[0] == '#') continue;
        int iR = iRz % nRND;
        int iz = iRz / nRND;
        if (RstND + iR*dRND == 1000*atof(words[0]) && zstND + iz*dzND == 1000*atof(words[1])) {
            logrhoNDs[iz][iR]    = log10(atof(words[2]));
            vphiNDs[iz][iR]      = atof(words[3]);
            logsigvNDs[iz][iR][0]= log10(atof(words[4]));
            logsigvNDs[iz][iR][1]= log10(atof(words[5]));
            logsigvNDs[iz][iR][2]= log10(atof(words[6]));
            corRzNDs[iz][iR]     = atof(words[7]);
        } else {
            printf("store_NSDmoments: index mismatch at iRz=%d\n", iRz);
        }
        iRz++;
    }
    fclose(fp);
}

// =========================================================
// IMF and stellar evolution
// =========================================================

void store_IMF_nBs(int B,
                   double *logMass, double *PlogM, double *PlogM_cum_norm, int *imptiles,
                   double M0, double M1, double M2, double M3, double Ml, double Mu,
                   double alpha1, double alpha2, double alpha3, double alpha4, double alpha0)
{
  // Verbatim from genulens.cpp store_IMF_nBs (Koshimoto, Baba & Bennett 2021)
  double *PlogM_cum, *Mass, *PMlogM_cum, *PMlogM_cum_norm;
  Mass            = (double*)calloc(nm+1, sizeof(double));
  PlogM_cum       = (double*)calloc(nm+1, sizeof(double));
  PMlogM_cum      = (double*)calloc(nm+1, sizeof(double));
  PMlogM_cum_norm = (double*)calloc(nm+1, sizeof(double));
  logMst = log10(Ml);
  dlogM  = (double)(log10(Mu)-logMst)/nm;
  for (int i = 0; i <= nm; i++) {
    double Mp  = i*dlogM + logMst;
    logMass[i] = Mp;
    Mass[i]    = pow(10, Mp);
    double alpha = (Mass[i] < M3) ? alpha4 : (Mass[i] < M2) ? alpha3
                 : (Mass[i] < M1) ? alpha2 : (Mass[i] < M0) ? alpha1 : alpha0;
    double temp00 = pow(M0, alpha0+1.), temp01 = pow(M0, alpha1+1.);
    double temp11 = pow(M1, alpha1+1.), temp12 = pow(M1, alpha2+1.);
    double temp22 = pow(M2, alpha2+1.), temp23 = pow(M2, alpha3+1.);
    double temp33 = pow(M3, alpha3+1.), temp34 = pow(M3, alpha4+1.);
    double dPlogM = 1;
    if (Mass[i] < M0) dPlogM = temp01/temp00;
    if (Mass[i] < M1) dPlogM = temp12/temp11*dPlogM;
    if (Mass[i] < M2) dPlogM = temp23/temp22*dPlogM;
    if (Mass[i] < M3) dPlogM = temp34/temp33*dPlogM;
    double templogMF = pow(Mass[i], alpha+1.);
    PlogM[i] = templogMF/dPlogM;
    if (i >= 1) {
      PlogM_cum[i]  = 0.5*(PlogM[i]+PlogM[i-1])*dlogM  + PlogM_cum[i-1];
      PMlogM_cum[i] = 0.5*(Mass[i]*PlogM[i]+Mass[i-1]*PlogM[i-1])*dlogM + PMlogM_cum[i-1];
    } else {
      PlogM_cum[i] = 0.0;
      PMlogM_cum[i] = 0.0;
    }
  }
  for (int i = 0; i <= nm; i++) {
    PlogM_cum_norm[i]  = PlogM_cum[i]  / PlogM_cum[nm];
    PMlogM_cum_norm[i] = PMlogM_cum[i] / PMlogM_cum[nm];
    PlogM[i] /= PlogM_cum[nm];
    int intp = (int)(PlogM_cum_norm[i]*20);
    if (imptiles[intp] == 0) imptiles[intp] = (intp == 0) ? 1 : (int)(i+0.5);
  }
  if (B == 0) { free(Mass); free(PlogM_cum); free(PMlogM_cum); free(PMlogM_cum_norm); return; }

  // ageMloss: cumulative average WD mass / initial mass, integrated from HIGH to LOW mass
  double *ageMloss = (double*)calloc(nm+1, sizeof(double));
  double cumMwt = 0, cumWDwt = 0;
  for (int i = nm; i >= 0; i--) {
    double pout[2] = {};
    double M  = pow(10, logMass[i]);
    double wt = PlogM[i];
    Mini2Mrem(pout, M, 1);
    double MWD = pout[0];
    cumMwt  += M   * wt;
    cumWDwt += MWD * wt;
    ageMloss[i] = cumWDwt / cumMwt;
  }

  // Read Minidie.dat
  char file1[] = "input_files/Minidie.dat";
  char line[1000]; char *words[100]; FILE *fp;
  double MRGstD[250], MRGenD[250], MRGstB[50], MRGenB[50], MRGstND[10], MRGenND[10];
  if ((fp = fopen(file1, "r")) == NULL) { printf("can't open %s\n", file1); exit(1); }
  nageD = 0; nageB = 0; nageND = 0;
  while (fgets(line, 1000, fp) != NULL) {
    split((char*)" ", line, words);
    if (*words[0] == '#') continue;
    if (*words[0] == 'N') {
      agesND[nageND]    = (int)atof(words[1]);
      MinidieND[nageND] = atof(words[2]);
      MRGstND[nageND]   = atof(words[3]);
      MRGenND[nageND]   = atof(words[4]);
      nageND++;
    } else if (*words[0] == 'B') {
      agesB[nageB]    = (int)atof(words[1]);
      MinidieB[nageB] = atof(words[2]);
      MRGstB[nageB]   = atof(words[3]);
      MRGenB[nageB]   = atof(words[4]);
      nageB++;
    } else {
      agesD[nageD]    = (int)atof(words[0]);
      MinidieD[nageD] = atof(words[1]);
      MRGstD[nageD]   = atof(words[2]);
      MRGenD[nageD]   = atof(words[3]);
      nageD++;
    }
  }
  fclose(fp);

  // Thin disk normalization
  double gamma = 1.0/tSFR;
  int agest = 1, ageen = 1000;
  int iages[7] = {15,100,200,300,500,700,1000};
  double wt_D[7]={}, wtWD_D[7]={}, sumM_D[7]={}, sumMWD_D[7]={};
  double sumstars_D[7]={}, sumWDs_D[7]={}, sumRGs_D[7]={};
  for (int i = agest; i <= ageen; i++) {
    int itmp = (int)((i - agesD[0])/(agesD[1] - agesD[0]) + 0.5);
    if (itmp < 0) itmp = 0;
    double logMdie = log10(MinidieD[itmp]);
    double logMRG1 = log10(MRGstD[itmp]);
    double logMRG2 = log10(MRGenD[itmp]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG  = PRG2 - PRG1;
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1-PM)*aveMloss;
    double PWD  = (1-P);
    double wtSFR = exp(-gamma*(ageen-i)*0.01);
    P    *= wtSFR; PWD *= wtSFR; PM   *= wtSFR; PMWD *= wtSFR; PRG  *= wtSFR;
    int idisk = (i <= iages[0]) ? 0 : (i <= iages[1]) ? 1 : (i <= iages[2]) ? 2
              : (i <= iages[3]) ? 3 : (i <= iages[4]) ? 4 : (i <= iages[5]) ? 5
              : (i <= iages[6]) ? 6 : 0;
    wt_D[idisk]      += PM;
    wtWD_D[idisk]    += PMWD;
    sumM_D[idisk]    += PM   * PMlogM_cum[nm];
    sumMWD_D[idisk]  += PMWD * PMlogM_cum[nm];
    sumstars_D[idisk]+= P    * PlogM_cum[nm];
    sumWDs_D[idisk]  += PWD  * PlogM_cum[nm];
    sumRGs_D[idisk]  += PRG  * PlogM_cum[nm];
  }
  // Compute rho0d, n0MSd, n0d, n0RGd for thin disk (i=0-6), thick disk (i=7), halo (i=8 in loop)
  double rho0thinMS = 0, rho0thinWD = 0;
  double Sig2rho[8]={}, aveMMS_D[8]={}, aveMWD_D[8]={}, nfracRG_D[8]={};
  double rhoT0 = rhot0 * 0.04;
  for (int i = 0; i < 8; i++) {
    Sig2rho[i] = 0.5/zd[i];
    if (i < 7) {
      int rd = (i == 0) ? Rd[0] : Rd[1];
      aveMMS_D[i]  = sumM_D[i]   / sumstars_D[i];
      aveMWD_D[i]  = sumMWD_D[i] / sumWDs_D[i];
      nfracRG_D[i] = sumRGs_D[i] / sumstars_D[i];
      rho0thinMS  += exp(-R0/rd) * wt_D[i]/rd * Sig2rho[i];
      rho0thinWD  += exp(-R0/rd) * wtWD_D[i]/rd * Sig2rho[i];
    }
  }
  for (int i = 0; i < 9; i++) {
    int rd = (i==0)?Rd[0]:(i<7)?Rd[1]:(i==7)?Rd[2]:0;
    if (i < 7) {
      double norm  = rhot0/rho0thinMS;
      double rhoMS = norm * exp(-R0/rd) * wt_D[i]/rd   * Sig2rho[i];
      double rhoWD = norm * exp(-R0/rd) * wtWD_D[i]/rd * Sig2rho[i];
      rho0d[i] = rhoMS + rhoWD;
      n0MSd[i] = rhoMS / aveMMS_D[i];
      double n0WD = rhoWD / aveMWD_D[i];
      n0d[i]   = n0MSd[i] + n0WD;
      n0RGd[i] = n0MSd[i] * nfracRG_D[i];
    } else if (i == 7) { // Thick disk
      double logMdie  = log10(MinidieD[nageD-2]);
      double logMRG1  = log10(MRGstD[nageD-2]);
      double logMRG2  = log10(MRGenD[nageD-2]);
      double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
      double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
      double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
      double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
      double PRG  = PRG2 - PRG1;
      double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
      double PMWD = (1-PM)*aveMloss;
      double PWD  = (1-P);
      double aveMMS = PM  *PMlogM_cum[nm]/P  /PlogM_cum[nm];
      double aveMWD = PMWD*PMlogM_cum[nm]/PWD/PlogM_cum[nm];
      double norm = rhoT0/PM;
      rho0d[7] = rhoT0 + norm*PMWD;
      n0MSd[7] = rhoT0/aveMMS;
      n0d[7]   = n0MSd[7] + norm*PMWD/aveMWD;
      n0RGd[7] = n0MSd[7]*PRG/P;
    } else { // Stellar halo
      double logMdie  = log10(MinidieD[nageD-1]);
      double logMRG1  = log10(MRGstD[nageD-1]);
      double logMRG2  = log10(MRGenD[nageD-1]);
      double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
      double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
      double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
      double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
      double PRG  = PRG2 - PRG1;
      double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
      double PMWD = (1-PM)*aveMloss;
      double PWD  = (1-P);
      double aveMMS = PM  *PMlogM_cum[nm]/P  /PlogM_cum[nm];
      double aveMWD = PMWD*PMlogM_cum[nm]/PWD/PlogM_cum[nm];
      double norm = rho0SHMS/PM;
      rho0SH  = rho0SHMS + norm*PMWD;
      n0MSSH  = rho0SHMS/aveMMS;
      n0SH    = n0MSSH + norm*PMWD/aveMWD;
      n0RGSH  = n0MSSH*PRG/P;
    }
  }
  // Bulge
  double wt_B=0,wtWD_B=0,sumM_B=0,sumMWD_B=0,sumstars_B=0,sumWDs_B=0,sumRGs_B=0;
  for (int i = 0; i < nageB; i++) {
    double tau    = 0.01*agesB[i];
    double wtSFR  = (tau-mageB)/sageB;
    wtSFR = exp(-0.5*wtSFR*wtSFR);
    double logMdie  = log10(MinidieB[i]);
    double logMRG1  = log10(MRGstB[i]);
    double logMRG2  = log10(MRGenB[i]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG  = PRG2 - PRG1;
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1-PM)*aveMloss;
    double PWD  = (1-P);
    P    *= wtSFR; PWD  *= wtSFR; PM   *= wtSFR; PMWD *= wtSFR; PRG  *= wtSFR;
    wt_B      += PM;  wtWD_B    += PMWD;
    sumM_B    += PM  *PMlogM_cum[nm]; sumMWD_B += PMWD*PMlogM_cum[nm];
    sumstars_B+= P   *PlogM_cum[nm];  sumWDs_B += PWD *PlogM_cum[nm];
    sumRGs_B  += PRG *PlogM_cum[nm];
  }
  double aveMMS_b = sumM_B   / sumstars_B;
  double aveMWD_b = sumMWD_B / sumWDs_B;
  m2nb_MS  = 1.0/aveMMS_b;
  m2nb_WD  = 1.0/aveMWD_b;
  nMS2nRGb = sumRGs_B/sumstars_B;
  fb_MS    = wt_B/(wt_B+wtWD_B);
  // NSD
  double wt_ND=0,wtWD_ND=0,sumM_ND=0,sumMWD_ND=0,sumstars_ND=0,sumWDs_ND=0,sumRGs_ND=0;
  for (int i = 0; i < nageND; i++) {
    double tau   = 0.01*agesND[i];
    double wtSFR = (tau-mageND)/sageND;
    wtSFR = exp(-0.5*wtSFR*wtSFR);
    double logMdie  = log10(MinidieND[i]);
    double logMRG1  = log10(MRGstND[i]);
    double logMRG2  = log10(MRGenND[i]);
    double PM   = interp_x(nm+1, PMlogM_cum_norm, logMst, dlogM, logMdie);
    double P    = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMdie);
    double PRG1 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG1);
    double PRG2 = interp_x(nm+1, PlogM_cum_norm,  logMst, dlogM, logMRG2);
    double PRG  = PRG2 - PRG1;
    double aveMloss = interp_x(nm+1, ageMloss, logMst, dlogM, logMdie);
    double PMWD = (1-PM)*aveMloss;
    double PWD  = (1-P);
    P    *= wtSFR; PWD  *= wtSFR; PM   *= wtSFR; PMWD *= wtSFR; PRG  *= wtSFR;
    wt_ND      += PM;  wtWD_ND    += PMWD;
    sumM_ND    += PM  *PMlogM_cum[nm]; sumMWD_ND += PMWD*PMlogM_cum[nm];
    sumstars_ND+= P   *PlogM_cum[nm];  sumWDs_ND += PWD *PlogM_cum[nm];
    sumRGs_ND  += PRG *PlogM_cum[nm];
  }
  double aveMMS_nd = sumM_ND   / sumstars_ND;
  double aveMWD_nd = sumMWD_ND / sumWDs_ND;
  m2nND_MS  = 1.0/aveMMS_nd;
  m2nND_WD  = 1.0/aveMWD_nd;
  nMS2nRGND = sumRGs_ND/sumstars_ND;
  fND_MS    = wt_ND/(wt_ND+wtWD_ND);
  free(Mass); free(PlogM_cum); free(PMlogM_cum); free(PMlogM_cum_norm); free(ageMloss);
}

void Mini2Mrem(double *pout, double Mini, int mean)
{
    double Mrem, fREM;
    double PNS = (Mini < MiniWDmax) ? 0
               : (Mini < 15.0)  ? 1
               : (Mini < 17.8)  ? 0.679
               : (Mini < 18.5)  ? 0.833
               : (Mini < 21.7)  ? 0.500
               : (Mini < 25.2)  ? 0
               : (Mini < 27.5)  ? 0.652
               : (Mini < 60.0)  ? 0
               : 0.4;
    if (Mini < MiniWDmax) {
        Mrem = 0.109*Mini + 0.394;
        fREM = 1;
    } else {
        double MNS;
        do {
            MNS = (Mini < 13.0) ? 2.24 + 0.508*(Mini-14.75)
                                        + 0.125*(Mini-14.75)*(Mini-14.75)
                                        + 0.011*(Mini-14.75)*(Mini-14.75)*(Mini-14.75)
                : (Mini < 15.0) ?  0.123 + 0.112*Mini
                : (Mini < 17.8) ?  0.996 + 0.0384*Mini
                : (Mini < 18.5) ? -0.020 + 0.10*Mini
                : (Mini < 21.7 && mean == 0) ? 1.60 + 0.158*gasdev()
                : (Mini < 21.7 && mean == 1) ? 1.60
                : (Mini < 27.5) ?  3232.29 - 409.429*(Mini-2.619)
                                           + 17.2867*(Mini-2.619)*(Mini-2.619)
                                           - 0.24315*(Mini-2.619)*(Mini-2.619)*(Mini-2.619)
                : (mean == 0) ? 1.78 + 0.02*gasdev() : 1.78;
        } while (PNS > 0 && (MNS < MNSMIN || MNS > MNSMAX));
        double Mcore = (Mini < 42.21) ? -2.049 + 0.4140*Mini
                     : 5.697 + 7.8598e8*pow(Mini, -4.858);
        double Mall  = 15.52 - 0.3294*(Mini-25.97)
                            - 0.02121*(Mini-25.97)*(Mini-25.97)
                            + 0.003120*(Mini-25.97)*(Mini-25.97)*(Mini-25.97);
        double fej = (Mini < 42.21) ? 0.9 : 1.0;
        double MBH = fej*Mcore + (1-fej)*Mall;
        if (mean == 1) {
            Mrem = PNS*MNS + (1-PNS)*MBH;
            fREM = PNS*2 + (1-PNS)*3;
        } else {
            double ran = ran1();
            Mrem = (ran < PNS) ? MNS : MBH;
            fREM = (ran < PNS) ? 2   : 3;
        }
    }
    pout[0] = Mrem;
    pout[1] = fREM;
}

// =========================================================
// Kinematics -- implemented in galactic_kinematics.cpp
// =========================================================

// =========================================================
// Model initializer
// =========================================================

void init_galactic_model(int argc, char **argv, int need_kinematics)
{
    // --- RNG ---
    long seed = getOptioni(argc, argv, "seed", 1, 12304357);
    gsl_rng_env_setup();
    T_rng = gsl_rng_default;
    r_rng = gsl_rng_alloc(T_rng);
    gsl_rng_set(r_rng, (unsigned long)seed);

    // --- IMF parameters ---
    double M0_B     = getOptiond(argc,argv,"M0",     1, 1.0);
    double M1_B     = getOptiond(argc,argv,"M1",     1, 0.859770466578045);
    double M2_B     = getOptiond(argc,argv,"M2",     1, 0.08);
    double M3_B     = getOptiond(argc,argv,"M3",     1, 0.01);
    double Ml       = getOptiond(argc,argv,"Ml",     1, 0.001);
    double Mu       = getOptiond(argc,argv,"Mu",     1, 120.0);
    double alpha1_B = getOptiond(argc,argv,"alpha1", 1, -2.32279457078378);
    double alpha2_B = getOptiond(argc,argv,"alpha2", 1, -1.13449983242887);
    double alpha3_B = getOptiond(argc,argv,"alpha3", 1, -0.175862190587576);
    double alpha0_B = getOptiond(argc,argv,"alpha0", 1, alpha1_B);
    double alpha4_B = getOptiond(argc,argv,"alpha4", 1, alpha3_B);

    // --- Disk flags ---
    DISK   = getOptiond(argc,argv,"DISK",   1, 2);
    rhot0  = getOptiond(argc,argv,"rhot0",  1, 0.042);
    hDISK  = getOptiond(argc,argv,"hDISK",  1, 0);
    addX   = getOptiond(argc,argv,"addX",   1, 5);
    model  = getOptiond(argc,argv,"model",  1, 5);
    B14disk= getOptiond(argc,argv,"B14disk",1, 0);
    int B14bar = getOptiond(argc,argv,"B14bar",1, 0);

    // --- Bar geometry ---
    R0     = getOptiond(argc,argv,"R0",     1, 8160);
    thetaD = getOptiond(argc,argv,"thetaD", 1, 27);
    frho0b = getOptiond(argc,argv,"frho0b", 1, 0.839014514507754);
    Rc     = getOptiond(argc,argv,"Rc",     1, 2631.78535429573);
    zb_c   = getOptiond(argc,argv,"zb_c",   1, 1e+6);
    if (model >= 4 && model <= 8) {
        x0_1 = getOptiond(argc,argv,"x0", 1, 930.623146993329);
        y0_1 = getOptiond(argc,argv,"y0", 1, 370.784386649364);
        z0_1 = getOptiond(argc,argv,"z0", 1, 239.547516030578);
        C1   = getOptiond(argc,argv,"C1", 1, 1.20011972384328);
        C2   = getOptiond(argc,argv,"C2", 1, 4.09326795684828);
        C3   = getOptiond(argc,argv,"C3", 1, 1.0);
    }
    if (addX >= 5) {
        x0_X = getOptiond(argc,argv,"x0_X", 1, 278.027059842233);
        y0_X = getOptiond(argc,argv,"y0_X", 1, 176.318528789193);
        z0_X = getOptiond(argc,argv,"z0_X", 1, 286.791941602401);
        C1_X = getOptiond(argc,argv,"C1_X", 1, 1.3087131258784);
        C2_X = getOptiond(argc,argv,"C2_X", 1, 2.21745322869032);
        b_zX = getOptiond(argc,argv,"b_zX", 1, 1.37774815817195);
        fX   = getOptiond(argc,argv,"fX",   1, 1.43975636704683);
        Rc_X = getOptiond(argc,argv,"Rc_X", 1, 1301.63829617294);
    }
    b_zY = getOptiond(argc,argv,"b_zY", 1, 0);

    // --- Bar kinematics (stored globally for get_vxyz_ran) ---
    Omega_p   = getOptiond(argc,argv,"Omega_p",  1, 47.4105844018699);
    model_vb  = getOptiond(argc,argv,"model_vb", 1, 5);
    x0_vb     = getOptiond(argc,argv,"x0_vb",    1, 858.106595717275);
    y0_vb     = getOptiond(argc,argv,"y0_vb",    1, 3217.04987721548);
    z0_vb     = getOptiond(argc,argv,"z0_vb",    1, 950.690583433628);
    C1_vb     = getOptiond(argc,argv,"C1_vb",    1, 4.25236641149869);
    C2_vb     = getOptiond(argc,argv,"C2_vb",    1, 1.02531652066343);
    C3_vb     = getOptiond(argc,argv,"C3_vb",    1, 1.0);
    sigx_vb   = getOptiond(argc,argv,"sigx_vb",  1, 151.854794853683);
    sigy_vb   = getOptiond(argc,argv,"sigy_vb",  1, 78.0278905748233);
    sigz_vb   = getOptiond(argc,argv,"sigz_vb",  1, 81.9641955092164);
    sigx_vb0  = getOptiond(argc,argv,"sigx_vb0", 1, 63.9939241108675);
    sigy_vb0  = getOptiond(argc,argv,"sigy_vb0", 1, 75.8180486866697);
    sigz_vb0  = getOptiond(argc,argv,"sigz_vb0", 1, 71.2336430487113);
    vx_str    = getOptiond(argc,argv,"vx_str",   1, 43.0364707040617);
    y0_str    = getOptiond(argc,argv,"y0_str",   1, 406.558313420815);
    model_vbz = getOptiond(argc,argv,"model_vbz",1, 5);
    x0_vbz    = getOptiond(argc,argv,"x0_vbz",   1, 558.430182718529);
    y0_vbz    = getOptiond(argc,argv,"y0_vbz",   1, 2003.21703656302);
    z0_vbz    = getOptiond(argc,argv,"z0_vbz",   1, 3823.20855045157);
    C1_vbz    = getOptiond(argc,argv,"C1_vbz",   1, 3.71001266000693);
    C2_vbz    = getOptiond(argc,argv,"C2_vbz",   1, 1.07455173734341);
    C3_vbz    = getOptiond(argc,argv,"C3_vbz",   1, 1.0);

    // --- Disk kinematics ---
    hsigUt  = getOptiond(argc,argv,"hsigUt", 1, 14300);
    hsigWt  = getOptiond(argc,argv,"hsigWt", 1, 5900);
    hsigUT  = getOptiond(argc,argv,"hsigUT", 1, 180000);
    hsigWT  = getOptiond(argc,argv,"hsigWT", 1, 9400);
    betaU   = getOptiond(argc,argv,"betaU",  1, 0.32);
    betaW   = getOptiond(argc,argv,"betaW",  1, 0.77);
    sigU10d = getOptiond(argc,argv,"sigU10d",1, 42.0);
    sigW10d = getOptiond(argc,argv,"sigW10d",1, 24.4);
    sigU0td = getOptiond(argc,argv,"sigU0td",1, 75.0);
    sigW0td = getOptiond(argc,argv,"sigW0td",1, 49.2);

    // --- Named models (overrides) ---
    int EXE_fg0 = getOptiond(argc,argv,"EXE_fg0",1,0);
    int E_fg0   = getOptiond(argc,argv,"E_fg0",  1,0);
    int G_fg0   = getOptiond(argc,argv,"G_fg0",  1,0);
    int GXG_fg0 = getOptiond(argc,argv,"GXG_fg0",1,0);
    if (EXE_fg0 == 1) {
        model=5; addX=5;
        M0_B=1.0; M1_B=0.859770466578045; M2_B=0.08; M3_B=0.01;
        alpha1_B=-2.32279457078378; alpha2_B=-1.13449983242887; alpha3_B=-0.175862190587576;
        alpha0_B=alpha1_B; alpha4_B=alpha3_B;
        R0=8160; thetaD=27; frho0b=0.839014514507754; Rc=2631.78535429573;
        x0_1=930.623146993329; y0_1=370.784386649364; z0_1=239.547516030578;
        C1=1.20011972384328; C2=4.09326795684828; C3=1;
        model_vb=5; model_vbz=5;
        Omega_p=47.4105844018699; vx_str=43.0364707040617; y0_str=406.558313420815;
        sigx_vb=151.854794853683; sigy_vb=78.0278905748233; sigz_vb=81.9641955092164;
        sigx_vb0=63.9939241108675; sigy_vb0=75.8180486866697; sigz_vb0=71.2336430487113;
        x0_vb=858.106595717275; y0_vb=3217.04987721548; z0_vb=950.690583433628;
        C1_vb=4.25236641149869; C2_vb=1.02531652066343;
        x0_vbz=558.430182718529; y0_vbz=2003.21703656302; z0_vbz=3823.20855045157;
        C1_vbz=3.71001266000693; C2_vbz=1.07455173734341;
        x0_X=278.027059842233; y0_X=176.318528789193; z0_X=286.791941602401;
        C1_X=1.3087131258784; C2_X=2.21745322869032;
        b_zX=1.37774815817195; fX=1.43975636704683; Rc_X=1301.63829617294;
    }
    if (B14bar == 1) {
        M0_B=1.0; M1_B=0.7; M2_B=0.08; M3_B=0.01;
        alpha1_B=-2.0; alpha2_B=-1.3; alpha3_B=-0.5;
        alpha0_B=alpha1_B; alpha4_B=alpha3_B;
        model=6; addX=0; B14vbar=1;
        R0=8200; thetaD=20; x0_1=1580; y0_1=620; z0_1=430; Rc=2400; C1=2; C2=4; C3=1;
        frho0b=1.173;
        Omega_p=50; sigx_vb=114; sigy_vb=103.8; sigz_vb=96.4;
        x0_vb=y0_vb=z0_vb=500000; x0_vbz=y0_vbz=z0_vbz=500000;
    }
    if (B14disk == 1) {
        vxsun=-12.7; vysun=218+24; vzsun=7.25; DISK=1;
    }

    // --- Derived geometry ---
    costheta = cos(thetaD/180.0*PI);
    sintheta = sin(thetaD/180.0*PI);

    // --- SgrA* position ---
    int CenSgrA = getOptioni(argc,argv,"CenSgrA",1,1);
    if (CenSgrA == 1) {
        double lSgrA = -0.0560, bSgrA = -0.0462; // deg
        Dlb2xyz(R0, lSgrA, bSgrA, R0, xyzSgrA);
    }

    // --- Stellar halo ---
    SH      = getOptiond(argc,argv,"SH",       1, 1);
    rho0SHMS= getOptiond(argc,argv,"rho0SHMS", 1, 9.32e-06);
    sigU_SH = getOptiond(argc,argv,"sigU_SH",  1, 131);
    sigV_SH = getOptiond(argc,argv,"sigV_SH",  1, 106);
    sigW_SH = getOptiond(argc,argv,"sigW_SH",  1, 85);
    if (SH == 0) rho0SHMS = 0;

    // --- Disk normalization at solar position ---
    y0d[0] = (DISK==1) ? exp(-R0/Rd[0] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[0]);
    y0d[1] = (DISK==1) ? exp(-R0/Rd[1] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[1]);
    y0d[2] = (DISK==1) ? exp(-R0/Rd[2] - pow((double)Rh/R0, nh)) : exp(-R0/Rd[2]);

    // --- IMF and normalization factors ---
    nm = 1000;
    g_logMass        = (double*)calloc(nm+1, sizeof(double));
    g_PlogM          = (double*)calloc(nm+1, sizeof(double));
    g_PlogM_cum_norm = (double*)calloc(nm+1, sizeof(double));
    int    *imptiles_B       = (int*)calloc(22, sizeof(int));
    store_IMF_nBs(1, g_logMass, g_PlogM, g_PlogM_cum_norm, imptiles_B,
                  M0_B, M1_B, M2_B, M3_B, Ml, Mu,
                  alpha1_B, alpha2_B, alpha3_B, alpha4_B, alpha0_B);
    free(imptiles_B);
    // g_logMass, g_PlogM, g_PlogM_cum_norm are kept for tools that need the IMF table

    // --- Bulge normalization (rho0b via crude_integrate) ---
    double massVVVbox = crude_integrate(2200, 1400, 1200, 15);
    double massentire = crude_integrate(6000, 3000, 3000, 30);
    double MVVVd = 0;
    for (int i = 0; i < 8; i++) {
        int rd = (i==0)?Rd[0]:(i<7)?Rd[1]:Rd[2];
        double zdtmp = (hDISK==0) ? zd[i] : zd45[i];
        double MVVVtmp = rho0d[i] * exp((R0-(double)Rdbreak)/rd) * 2200*2 * 1400*2 * zd[i]/zdtmp;
        double ztmp = 1200.0/zdtmp;
        MVVVtmp *= (i<7) ? 2*zdtmp*(exp(2*ztmp)-1)/(exp(2*ztmp)+1)
                         : 2*zdtmp*(1 - exp(-ztmp));
        MVVVd += MVVVtmp;
    }
    double MVVVP17 = 1.32e+10;
    rho0b  = (frho0b*MVVVP17 - MVVVd) / massVVVbox;
    n0MSb  = rho0b * fb_MS * m2nb_MS;
    n0RGb  = n0MSb * nMS2nRGb;
    n0b    = n0MSb + rho0b * (1-fb_MS) * m2nb_WD;

    // --- NSD setup ---
    double lSIMU = getOptiond(argc,argv,"l", 1, 1.0);
    double bSIMU = getOptiond(argc,argv,"b", 1, -3.9);
    if (fabs(lSIMU) < 5 && fabs(bSIMU) < 2) ND = 3;
    ND = (int)getOptiond(argc,argv,"NSD", 1, ND);
    if (ND > 0) ND = 3;
    double MND = 1e9; // default (unused when ND==3)
    if (ND == 1) { MND=2e9; x0ND=250; y0ND=125; z0ND=50; }
    if (ND == 2) { MND=7e8; x0ND=74;  y0ND=74;  z0ND=26; }
    x0ND = (int)getOptiond(argc,argv,"x0ND",1,x0ND);
    y0ND = (int)getOptiond(argc,argv,"y0ND",1,y0ND);
    z0ND = (int)getOptiond(argc,argv,"z0ND",1,z0ND);
    MND  = getOptiond(argc,argv,"MND", 1, MND);
    if (ND > 0) {
        rho0ND = (ND==3) ? 1 : 0.25*MND/PI/x0ND/y0ND/z0ND;
        n0MSND = rho0ND * fND_MS * m2nND_MS;
        n0RGND = n0MSND * nMS2nRGND;
        n0ND   = n0MSND + rho0ND*(1-fND_MS)*m2nND_WD;
    }
    nzND = (int)((zenND - zstND)/dzND + 1.5);
    nRND = (int)((RenND - RstND)/dRND + 1.5);
    if (ND == 3) {
        logrhoNDs  = (double**)malloc(sizeof(double*)*nzND);
        vphiNDs    = (double**)malloc(sizeof(double*)*nzND);
        corRzNDs   = (double**)malloc(sizeof(double*)*nzND);
        logsigvNDs = (double***)malloc(sizeof(double**)*nzND);
        for (int i = 0; i < nzND; i++) {
            logrhoNDs[i]  = (double*)calloc(nRND, sizeof(double));
            vphiNDs[i]    = (double*)calloc(nRND, sizeof(double));
            corRzNDs[i]   = (double*)calloc(nRND, sizeof(double));
            logsigvNDs[i] = (double**)malloc(sizeof(double*)*nRND);
            for (int j = 0; j < nRND; j++) {
                logsigvNDs[i][j] = (double*)calloc(3, sizeof(double));
            }
        }
        char *fileND = (char*)"input_files/NSD_moments.dat";
        store_NSDmoments(fileND);
    }

    // --- LOS coordinate for idata=0 ---
    lDs = (double*)malloc(sizeof(double));
    bDs = (double*)malloc(sizeof(double));
    lDs[0] = lSIMU;
    bDs[0] = bSIMU;

    // --- Kinematics (Phase 3 only) ---
    if (need_kinematics) {
        char *fileVc = (char*)"input_files/Rotcurve_BG16.dat";
        store_cumuP_Shu(fileVc);
    }
}
