/* galactic_kinematics.cpp
 * Shu distribution-function velocity sampling for pre_gapmoe tools.
 * Implements store_cumuP_Shu() and get_vxyz_ran() declared in galactic_model.h.
 * All functions are verbatim or lightly adapted from genulens.cpp
 * (Koshimoto, Baba & Bennett 2021).  NR functions replaced with GSL equivalents
 * via the ran1()/gasdev() wrappers already in galactic_model.cpp.
 */

#include "galactic_model.h"
#include "option.h"

#define vescd 550.0   // escape velocity disk  [km/s]
#define vescb 600.0   // escape velocity bulge [km/s]

// =========================================================
// Local helpers (not exposed in galactic_model.h)
// =========================================================

static double getcumu2xist(int n, double *x, double *F, double *f,
                            double Freq, int ist, int inv)
{
    double Fmax = F[n-1], Fmin = F[0];
    if (Fmin > Freq) return 0;
    if (Fmax < Freq) return 0;
    if (ist < 1) ist = 1;
    if (inv == 0) {
        for (int i = ist; i < n; i++) {
            if ((F[i] <= Freq && F[i-1] > Freq) || (F[i] >= Freq && F[i-1] < Freq)) {
                double a = 0.5*(f[i]-f[i-1])/(x[i]-x[i-1]);
                double b = f[i-1] - 2*a*x[i-1];
                double c = a*x[i-1]*x[i-1] - f[i-1]*x[i-1] + F[i-1] - Freq;
                return (a != 0) ? (-b + sqrt(b*b - 4*a*c)) * 0.5/a
                                : (x[i]-x[i-1])/(F[i]-F[i-1])*(Freq-F[i-1]) + x[i-1];
            }
        }
    } else {
        for (int i = ist; i > 0; i--) {
            if ((F[i] <= Freq && F[i-1] > Freq) || (F[i] >= Freq && F[i-1] < Freq)) {
                double a = 0.5*(f[i]-f[i-1])/(x[i]-x[i-1]);
                double b = f[i-1] - 2*a*x[i-1];
                double c = a*x[i-1]*x[i-1] - f[i-1]*x[i-1] + F[i-1] - Freq;
                return (a != 0) ? (-b + sqrt(b*b - 4*a*c)) * 0.5/a
                                : (x[i]-x[i-1])/(F[i]-F[i-1])*(Freq-F[i-1]) + x[i-1];
            }
        }
    }
    return 0;
}

static double calc_gc(double c)
{
    if (c < 0.5) return 0;
    if (c < 10) {
        double c2 = c - 0.5;
        return exp(c) * tgamma(c2) / (2 * pow(c, c2));
    }
    return sqrt(0.5*PI / (c - 0.913));
}

static double calc_SigRg(double Rg, double hsigU, int rd, double a0)
{
    double k=31.53, a=0.6719, b=0.2743;
    double c1=3.822, c2=0.524, c3=0.00567, c4=2.13;
    double q  = (double)rd/hsigU;
    double Rgmax = c1*rd/(1+q/c2);
    double x  = Rg/Rgmax;
    double s  = k*exp(-x/b)*((x/a)*(x/a) - 1);
    return 0.5*exp(-Rg/rd)/PI - c3*pow(a0,c4)*s;
}

static double calc_faca(double Rg, double hsigU, int rd, double a0)
{
    double q = (double)rd/hsigU;
    double bunsi = 0.25*pow(a0, 2.04);
    double bumbo = pow(q, 0.49);
    double as[12] = {-0.028476,-1.4518,12.492,-21.842,19.130,-10.175,3.5214,
                     -0.81052,0.12311,-0.011851,0.00065476,-1.5809e-05};
    double x = Rg*q/rd;
    double fp = as[0];
    double xp = 1;
    for (int j = 1; j < 12; j++) { xp *= x; fp += as[j]*xp; }
    return 1 - bunsi/bumbo * fp;
}

static double calc_PRRg(int R, int z, double fg, double sigU0, double hsigU, int rd)
{
    if (fg <= 0) return 0;
    double Rg  = R*fg;
    double vc  = getx2y(nVcs, Rcs, Vcs, Rg) / (1 + 0.0374*pow(0.001*abs(z), 1.34));
    double a0  = sigU0/vc * exp(R0/hsigU);
    double a   = sigU0/vc * exp(-(Rg - R0)/hsigU) * calc_faca(Rg, hsigU, rd, a0);
    double c   = 0.5/a/a;
    if (c <= 0.5) return 0;
    double SigRg = calc_SigRg(Rg, hsigU, rd, a0);
    double gc    = calc_gc(c);
    double x     = c*(2*log(fg) + 1 - fg*fg);
    double PRRg  = SigRg * exp(x)/gc;
    return (PRRg < 0) ? 0 : PRRg;
}

static void calc_dpdfg(double *pout, int R, int z, double fg1, double sigU0, double hsigU, int rd)
{
    double dfg  = 0.001;
    double fg2  = fg1 + dfg;
    double P1   = calc_PRRg(R, z, fg1, sigU0, hsigU, rd);
    double P2   = calc_PRRg(R, z, fg2, sigU0, hsigU, rd);
    double dP   = (P2 <= 0 || P1 <= 0) ? 0 : (P2-P1)/dfg;
    if (P2 <= 0 || P1 <= 0) P1 = 0;
    pout[0] = dP;
    pout[1] = P1;
}

static void get_PRRGmax2(double *pout, int R, int z, double fg1,
                          double sigU0, double hsigU, int rd)
{
    if (fg1 < 1) fg1 = 1;
    double dfg = 0.001;
    double fgc = 1e3, Pmax = 1e-200, dPdfgc = 0;
    double Ptmp = 0;
    double fg, fg2, fg3, fg4, dPdfg1, dPdfg2, d2Pdfg, dPdfg3, dPdfg4, d2Pdfg2, P1, P2, P3, P4;
    double jj;
    int nj = 0, ntry = 0, sw = 0;
    double pout1[2]={}, pout2[2]={};
    if (hsigU/rd/sigU0 < 0.1) {
        for (fg = 0.15; fg < 1.0; fg += 0.05) {
            P1 = calc_PRRg(R, z, fg, sigU0, hsigU, rd);
            if (P1 > Ptmp) { Ptmp = P1; fg1 = fg; }
        }
    }
    while (1) {
        int ncalc = 0;
        for (int j = 0; j < 3; j++) {
            fg2 = fg1 + dfg;
            calc_dpdfg(pout1, R, z, fg1, sigU0, hsigU, rd);
            calc_dpdfg(pout2, R, z, fg2, sigU0, hsigU, rd);
            dPdfg1 = pout1[0]; dPdfg2 = pout2[0];
            P1 = pout1[1];     P2 = pout2[1];
            d2Pdfg = (dPdfg2 - dPdfg1)/dfg;
            if (P1 > Pmax) { fgc = fg1; dPdfgc = dPdfg1; Pmax = P1; }
            ncalc++;
            if (ncalc > 15) {
                if (nj > 0) break;
                else if (ntry < 2) {
                    if (fgc > 900) fgc = (ntry == 0) ? fg1 : 0.9;
                    fg1 = (ntry == 0) ? fgc - 0.4 : fgc + 0.4;
                    if (fg1 < 0) fg1 = 0.2*ran1();
                    ncalc = 0; j = -1; ntry++; continue;
                } else break;
            }
            if (j == 2 && fabs(dPdfgc/Pmax) > 0.1) {
                nj++;
                fg1 = (dPdfgc > 0) ? fgc + 0.05/nj*ran1() : fgc - 0.05/nj*ran1();
                j = -1; continue;
            }
            if (dPdfg1 == 0) {
                jj = (dPdfgc == 0) ? 0.5 : 0.2*ran1();
                fg1 = (fg1 < fgc) ? fg1 + jj : fg1 - jj;
                j = -1; continue;
            }
            if (d2Pdfg > 0 && dPdfg1 < 0) {
                fg3 = fg2 + 0.04; fg4 = fg3 + dfg;
                calc_dpdfg(pout1, R, z, fg3, sigU0, hsigU, rd);
                calc_dpdfg(pout2, R, z, fg4, sigU0, hsigU, rd);
                dPdfg3 = pout1[0]; dPdfg4 = pout2[0];
                d2Pdfg2 = (dPdfg4-dPdfg3)/dfg;
                if (d2Pdfg2 > 0 || dPdfg3 == 0) { fg1 -= (0.02 + 0.10*ran1()); j = -1; continue; }
                else d2Pdfg = d2Pdfg2;
            }
            if (d2Pdfg > 0 && dPdfg1 > 0) {
                fg3 = fg1 - 0.04; fg4 = fg3 + dfg;
                calc_dpdfg(pout1, R, z, fg3, sigU0, hsigU, rd);
                calc_dpdfg(pout2, R, z, fg4, sigU0, hsigU, rd);
                dPdfg3 = pout1[0]; dPdfg4 = pout2[0];
                d2Pdfg2 = (dPdfg4-dPdfg3)/dfg;
                if (d2Pdfg2 > 0 || dPdfg3 == 0) { fg1 += (0.02 + 0.10*ran1()); j = -1; continue; }
                else d2Pdfg = d2Pdfg2;
            }
            if (d2Pdfg != 0) fg1 = fg1 - dPdfg1/d2Pdfg;
            if (fg1 < 0) fg1 = 0.1;
            if (fabs(dPdfg1/d2Pdfg) > 0.5) {
                jj = (dPdfgc > 0) ? 0.10 : -0.10;
                fg1 = fgc + jj*ran1();
                j = -1; continue;
            }
        }
        ncalc = 0; sw = 0;
        for (fg1 = fgc-0.2; fg1 > 0.1; fg1 -= 0.2) {
            P1 = calc_PRRg(R, z, fg1, sigU0, hsigU, rd);
            ncalc++;
            if (P1 > Pmax*1.05) { Pmax = P1; fgc = fg1; sw = 1; }
            if (P1/Pmax < 1e-02) break;
        }
        if (sw == 1) { fg1 = fgc; continue; }
        sw = 0;
        for (fg2 = fgc+0.2; fg2 < 4.0; fg2 += 0.2) {
            P2 = calc_PRRg(R, z, fg2, sigU0, hsigU, rd);
            ncalc++;
            if (P2 > Pmax*1.05) { Pmax = P2; fgc = fg2; sw = 1; break; }
            if (P2/Pmax < 1e-02) break;
        }
        if (sw == 1) { fg1 = fgc; continue; }
        if (fg1 < 0) fg1 = 0.1;
        pout[0] = Pmax; pout[1] = fg1; pout[2] = fg2; pout[3] = fgc;
        break;
    }
}

static void calc_sigvb(double xb, double yb, double zb, double *sigvbs)
{
    double xn, yn, zn, Rs, rs, facsig, facsigz = 0;
    xn = fabs(xb/x0_vb); yn = fabs(yb/y0_vb); zn = fabs(zb/z0_vb);
    Rs = pow(pow(xn,C1_vb) + pow(yn,C1_vb), 1.0/C1_vb);
    rs = pow(pow(Rs,C2_vb) + pow(zn,C2_vb), 1.0/C2_vb);
    if (rs == 0 && model_vb == 8) rs = 0.0001;
    facsig = (model_vb == 5) ? exp(-rs)
           : (model_vb == 6) ? exp(-0.5*rs*rs)
           : (model_vb == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2)
           : (model_vb == 4) ? exp(-pow(rs, C3_vb))
           : 0;
    if (model_vbz >= 4) {
        xn = fabs(xb/x0_vbz); yn = fabs(yb/y0_vbz); zn = fabs(zb/z0_vbz);
        Rs = pow(pow(xn,C1_vbz) + pow(yn,C1_vbz), 1.0/C1_vbz);
        rs = pow(pow(Rs,C2_vbz) + pow(zn,C2_vbz), 1.0/C2_vbz);
        if (rs == 0 && model_vbz == 8) rs = 0.0001;
        facsigz = (model_vbz == 5) ? exp(-rs)
                : (model_vbz == 6) ? exp(-0.5*rs*rs)
                : (model_vbz == 7) ? pow(2.0/(exp(rs)+exp(-rs)), 2)
                : (model_vbz == 4) ? exp(-pow(rs, C3_vbz))
                : 0;
    } else {
        facsigz = facsig;
    }
    sigvbs[0] = sigx_vb * facsig  + sigx_vb0;
    sigvbs[1] = sigy_vb * facsig  + sigy_vb0;
    sigvbs[2] = sigz_vb * facsigz + sigz_vb0;
}

// =========================================================
// Public API (declared in galactic_model.h)
// =========================================================

void store_cumuP_Shu(char *infile)
{
    // Read rotation curve
    FILE *fp;
    char line[1000];
    char *words[100];
    if (nVcs == 0) {
        if ((fp = fopen(infile, "r")) == NULL) { printf("can't open %s\n", infile); exit(1); }
        nVcs = 0;
        while (fgets(line, 1000, fp) != NULL) {
            split((char*)" ", line, words);
            if (*words[0] == '#') continue;
            Rcs[nVcs] = 1000*atof(words[0]);
            Vcs[nVcs] =      atof(words[1]);
            nVcs++;
        }
        fclose(fp);
    }

    // Allocate Shu DF table arrays
    int nfg   = 100;
    int nz    = (zenShu - zstShu)/dzShu + 1;
    int nR    = (RenShu - RstShu)/dRShu + 1;
    int ndisk = 8;
    fgsShu     = (double****)malloc(sizeof(double***) * nz);
    PRRgShus   = (double****)malloc(sizeof(double***) * nz);
    cumu_PRRgs = (double****)malloc(sizeof(double***) * nz);
    n_fgsShu   = (int***)malloc(sizeof(int**) * nz);
    kptiles    = (int****)malloc(sizeof(int***) * nz);
    for (int i = 0; i < nz; i++) {
        fgsShu[i]     = (double***)malloc(sizeof(double**) * nR);
        PRRgShus[i]   = (double***)malloc(sizeof(double**) * nR);
        cumu_PRRgs[i] = (double***)malloc(sizeof(double**) * nR);
        n_fgsShu[i]   = (int**)malloc(sizeof(int*) * nR);
        kptiles[i]    = (int***)malloc(sizeof(int**) * nR);
        for (int j = 0; j < nR; j++) {
            fgsShu[i][j]     = (double**)malloc(sizeof(double*) * ndisk);
            PRRgShus[i][j]   = (double**)malloc(sizeof(double*) * ndisk);
            cumu_PRRgs[i][j] = (double**)malloc(sizeof(double*) * ndisk);
            kptiles[i][j]    = (int**)malloc(sizeof(int*) * ndisk);
            n_fgsShu[i][j]   = (int*)calloc(ndisk, sizeof(int));
            for (int k = 0; k < ndisk; k++) {
                fgsShu[i][j][k]     = (double*)calloc(nfg, sizeof(double));
                PRRgShus[i][j][k]   = (double*)calloc(nfg, sizeof(double));
                cumu_PRRgs[i][j][k] = (double*)calloc(nfg, sizeof(double));
                kptiles[i][j][k]    = (int*)calloc(22, sizeof(int));
            }
        }
    }

    // Fill cumulative Shu DF table
    for (int z = zstShu; z <= zenShu; z += dzShu) {
        int iz = (z - zstShu)/dzShu;
        double facVcz = 1 + 0.0374*pow(0.001*abs(z), 1.34);
        for (int R = RstShu; R <= RenShu; R += dRShu) {
            int iR = (R - RstShu)/dRShu;
            double vcR = getx2y(nVcs, Rcs, Vcs, R);
            for (int idisk = 0; idisk < 8; idisk++) {
                double tau   = medtauds[idisk];
                double hsigU = (idisk < 7) ? hsigUt : hsigUT;
                int    rd    = (idisk == 0) ? Rd[0] : (idisk < 7) ? Rd[1] : Rd[2];
                double sigU0 = (idisk < 7) ? sigU10d*pow((tau+0.01)/10.01, betaU) : sigU0td;
                double Rgmin = R0 - hsigU*log(vcR/sigU0);
                if (Rgmin > R) Rgmin = R0 - hsigU*log(240.0/sigU0);
                double fgmin0 = Rgmin/R;
                double fg1 = (fgmin0 > 1.5) ? fgmin0 : 1;
                double pout[4] = {};
                get_PRRGmax2(pout, R, z, fg1, sigU0, hsigU, rd);
                double Pmax = pout[0], fgmin = pout[1], fgmax = pout[2], fgc = pout[3];
                if ((fgmin > 1 && R > 1000) || Pmax == 0)
                    printf("# PERROR!! get_PRRGmax2(pout, %5d, %4d, %.3f, %.2f, %.2f, %d)\n",
                           R, z, fg1, sigU0, hsigU, rd);
                int swerror = ((fgmin > 1 && R > 1000) || Pmax == 0) ? 1 : 0;
                double fg  = fgmin;
                double dfg0 = (fgc - fgmin)*0.025;
                int ifg = 0;
                double dfg = 0;
                while (fg <= fgmax) {
                    fgsShu[iz][iR][idisk][ifg]     = fg;
                    double PRRg = calc_PRRg(R, z, fg, sigU0, hsigU, rd);
                    PRRgShus[iz][iR][idisk][ifg]   = PRRg;
                    cumu_PRRgs[iz][iR][idisk][ifg] = (ifg == 0) ? 0
                        : cumu_PRRgs[iz][iR][idisk][ifg-1]
                          + 0.5*(PRRgShus[iz][iR][idisk][ifg-1] + PRRgShus[iz][iR][idisk][ifg])*dfg;
                    dfg = (PRRg/Pmax < 0.05) ? 4*dfg0
                        : (PRRg/Pmax < 0.25 || PRRg/Pmax > 0.7) ? dfg0
                        : 2*dfg0;
                    ifg++;
                    fg += dfg;
                }
                n_fgsShu[iz][iR][idisk] = ifg;
                double norm = cumu_PRRgs[iz][iR][idisk][ifg-1];
                for (int k = 0; k < ifg; k++) {
                    PRRgShus[iz][iR][idisk][k]   /= norm;
                    cumu_PRRgs[iz][iR][idisk][k] /= norm;
                    int intp = (int)(cumu_PRRgs[iz][iR][idisk][k]*20);
                    if (kptiles[iz][iR][idisk][intp] == 0)
                        kptiles[iz][iR][idisk][intp] = (intp == 0) ? 1 : (int)(k+0.5);
                }
                if (swerror == 1)
                    printf("# i=%d tau=%5.2f fg=%7.4f-%7.4f fgc=%6.4f Pmax=%.3e\n",
                           idisk, tau, fgmin, fgmax, fgc, Pmax);
            }
        }
    }
}

void get_vxyz_ran(double *vxyz, int i, double tau, double D, double lD, double bD)
{
    double xyz[3] = {};
    Dlb2xyz(D, lD, bD, R0, xyz);
    double x = xyz[0], y = xyz[1], z = xyz[2];
    double R = sqrt(x*x + y*y);
    double vx = 0, vy = 0, vz = 0;

    if (i < 8 && B14disk == 1) {
        double aveV = (i < 7) ? 218.0 : 170.0;
        double sigU = (i < 7) ?  39.9 :  67.0;
        double sigV = (i < 7) ?  27.9 :  51.0;
        double sigW = (i < 7) ?  19.1 :  42.0;
        do {
            double vR   =    0 + gasdev()*sigU;
            double vphi = aveV + gasdev()*sigV;
            vx = -vphi*y/R + vR*x/R;
            vy =  vphi*x/R + vR*y/R;
            vz =    0 + gasdev()*sigW;
        } while (vx*vx + vy*vy + vz*vz > vescd*vescd);

    } else if (i < 8) {
        double sigW0 = (i < 7) ? sigW10d*pow((tau+0.01)/10.01, betaW) : sigW0td;
        double sigU0 = (i < 7) ? sigU10d*pow((tau+0.01)/10.01, betaU) : sigU0td;
        double hsigW = (i < 7) ? hsigWt : hsigWT;
        double hsigU = (i < 7) ? hsigUt : hsigUT;
        double sigW  = sigW0*exp(-(R - R0)/hsigW);
        double sigU  = sigU0*exp(-(R - R0)/hsigU);
        int iz = (int)((fabs(z) - zstShu)/dzShu);
        int iR = (R > RstShu) ? (int)((R - RstShu)/dRShu) : 0;
        do {
            double ran = ran1();
            int inttmp = (int)(ran*20);
            int kst1 = 1;
            for (int itmp = inttmp; itmp > 0; itmp--)
                if (kst1 == 1) { kst1 = kptiles[iz][iR][i][itmp]; break; }
            double fg = getcumu2xist(n_fgsShu[iz][iR][i], fgsShu[iz][iR][i],
                                     cumu_PRRgs[iz][iR][i], PRRgShus[iz][iR][i],
                                     ran, kst1, 0);
            double Rg   = fg*R;
            double vc   = getx2y(nVcs, Rcs, Vcs, Rg) / (1 + 0.0374*pow(0.001*fabs(z), 1.34));
            double vphi = vc*fg;
            double vR   = gasdev()*sigU;
            vx = -vphi*y/R + vR*x/R;
            vy =  vphi*x/R + vR*y/R;
            vz = gasdev()*sigW;
        } while (vx*vx + vy*vy + vz*vz > vescd*vescd);

    } else if (i == 9 && ND == 3) {
        if (R > RenND || fabs(z) > zenND) {
            printf("ERROR: NSD comp exists where it must not. (R,z)=(%.0f,%.0f)\n", R, z);
            exit(1);
        }
        double as[4] = {};
        interp_xy_coeff(nzND, nRND, as, zstND, RstND, dzND, dRND, fabs(z), R);
        int iz0 = (int)((fabs(z) - zstND)/dzND);
        int iR0 = (int)((R - RstND)/dRND);
        double m_vphi=0, logsigphi=0, logsigR=0, logsigz=0, corRz=0;
        for (int j = 0; j < 4; j++) {
            int izj = (j == 0 || j == 2) ? iz0 : iz0+1;
            int iRj = (j == 0 || j == 1) ? iR0 : iR0+1;
            if (as[j] > 0) {
                m_vphi    += as[j]*vphiNDs[izj][iRj];
                logsigphi += as[j]*logsigvNDs[izj][iRj][0];
                logsigR   += as[j]*logsigvNDs[izj][iRj][1];
                logsigz   += as[j]*logsigvNDs[izj][iRj][2];
                corRz     += as[j]*corRzNDs[izj][iRj];
            }
        }
        double sigphi  = pow(10.0, logsigphi);
        double sigR    = pow(10.0, logsigR);
        double sigzv   = pow(10.0, logsigz);
        double facR    = sigzv/sigR * corRz;
        double sigz_R  = sigzv*sqrt(1 - corRz*corRz);
        do {
            double vphi = m_vphi + gasdev()*sigphi;
            double vR   = gasdev()*sigR;
            vx = -vphi*y/R + vR*x/R;
            vy =  vphi*x/R + vR*y/R;
            vz = facR*vR + gasdev()*sigz_R;
        } while (vx*vx + vy*vy + vz*vz > vescb*vescb);

    } else if (i == 10) {
        double vR   = sigU_SH * gasdev();
        double vphi = sigV_SH * gasdev();
        vx = -vphi*y/R + vR*x/R;
        vy =  vphi*x/R + vR*y/R;
        vz = sigW_SH * gasdev();

    } else { // bar & NSD (ND <= 2)
        double vrot = 0.001 * Omega_p * R;
        if (B14vbar == 1 && vrot > 90) vrot = 90;
        double xb =  x*costheta + y*sintheta;
        double yb = -x*sintheta + y*costheta;
        double zb =  z;
        double sigvbs[3] = {};
        double sigx, sigy, sigz;
        if (B14vbar == 1) {
            double Mvb_xrot = -vrot*y/R;
            double Mvb_yrot =  vrot*x/R;
            sigx = sqrt(sigx_vb*sigx_vb - Mvb_xrot*Mvb_xrot);
            sigy = sqrt(sigy_vb*sigy_vb - Mvb_yrot*Mvb_yrot);
            sigz = sigz_vb;
        } else {
            calc_sigvb(xb, yb, zb, sigvbs);
            sigx = sqrt(sigvbs[0]*sigvbs[0]*costheta*costheta + sigvbs[1]*sigvbs[1]*sintheta*sintheta);
            sigy = sqrt(sigvbs[0]*sigvbs[0]*sintheta*sintheta + sigvbs[1]*sigvbs[1]*costheta*costheta);
            sigz = sigvbs[2];
        }
        double avevxb = (yb > 0) ? -vx_str : vx_str;
        if (y0_str > 0) {
            double tmpyn = fabs(yb/y0_str);
            avevxb *= (1 - exp(-tmpyn*tmpyn));
        }
        do {
            vx = -vrot*y/R + avevxb*costheta + sigx*gasdev();
            vy =  vrot*x/R + avevxb*sintheta + sigy*gasdev();
            vz =                               sigz*gasdev();
        } while (vx*vx + vy*vy + vz*vz > vescb*vescb);
    }

    vxyz[0] = vx;
    vxyz[1] = vy;
    vxyz[2] = vz;
}
