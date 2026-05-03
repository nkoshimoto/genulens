#include "genulens/simulation/los_density_grid.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"

namespace genulens {

static const double a2toSig25BHs[5] = {0,  0.00279,  0.00985,  0.02054,  0.01901};
static const double a1toSig25BHs[5] = {0, -0.02023, -0.07548, -0.16585, -0.15082};
static const double a0toSig25BHs[5] = {1,  1.01342,  1.05382,  1.05898,  0.69889};

void LineOfSightDensityGrid::build(RunContext &ctx,
                                    const LineOfSightDensityGridConfig &cfg) {
    active_state = &ctx;

    ncomp_ = ncomp;
    const int Dmax = cfg.Dmax;

    nbin_ = (ND > 0 && fabs(cfg.l) < 0.05 && fabs(cfg.b) < 0.05) ? (int)(0.20*Dmax + 0.5)
          : (ND > 0 && fabs(cfg.l) < 0.10 && fabs(cfg.b) < 0.10) ? (int)(0.10*Dmax + 0.5)
          : (ND > 0)                                               ? (int)(0.04*Dmax + 0.5)
          :                                                           (int)(0.01*Dmax + 0.5);
    dD_    = (double)Dmax / nbin_;
    nallS_ = 0.0;

    D_.assign(nbin_ + 1, 0.0);
    cumu_rho_all_S_.assign(nbin_ + 1, 0.0);
    cumu_rho_all_L_.assign(nbin_ + 1, 0.0);
    rhoD_S_.assign(ncomp_, std::vector<double>(nbin_ + 1, 0.0));
    rhoD_L_.assign(ncomp_, std::vector<double>(nbin_ + 1, 0.0));
    cumu_rho_S_.assign(ncomp_, std::vector<double>(nbin_ + 1, 0.0));
    cumu_rho_L_.assign(ncomp_, std::vector<double>(nbin_ + 1, 0.0));
    fBH_.assign(ncomp_, std::vector<double>(nbin_ + 1, 1.0));
    ibinptiles_S_.assign(ncomp_, std::vector<int>(22, 0));
    ibinptiles_L_.assign(ncomp_, std::vector<int>(22, 0));

    void calc_rho_each(double D, int idata, double *rhos, double *xyz, double *xyb);
    double fLF_detect(double extI, double Imin, double Imax, int idisk);
    double fIVI_detect(double extI, double Imin, double Imax,
                       double extVI, double VImin, double VImax, int idisk);

    double rhos[11] = {}, xyz[3] = {}, xyb[2] = {};
    const int idata = 0;

    for (int ibin = 0; ibin <= nbin_; ibin++) {
        D_[ibin] = (double)ibin / nbin_ * Dmax;
        calc_rho_each(D_[ibin], idata, rhos, xyz, xyb);
        double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
        double z = xyz[2];

        if (ibin % cfg.npri == 0)
            printf("# %5.0f %5.0f %5.0f ", D_[ibin], R, z);

        double rhosum = 0;
        const auto dust = cfg.extinction->at_distance(D_[ibin]);
        double extI  = dust.i_band + cfg.extinction->distance_modulus_term(D_[ibin] + 0.1);
        double extVI = dust.color_vi;

        for (int i = 0; i < ncomp_; i++) {
            double fBHtmp = 1.0;

            if (cfg.UseSigBH == 1 && i < 9 && cfg.vkickBH > 25) {
                int ivkick = (int)floor(log(cfg.vkickBH / 25) / log(2) + 0.5);
                double corSigFac = a2toSig25BHs[ivkick]*R*R*1e-06
                                 + a1toSig25BHs[ivkick]*R*1e-03
                                 + a0toSig25BHs[ivkick];
                fBHtmp *= corSigFac;
            }

            int ien = (cfg.BHhb == 1) ? 9 : 8;
            if (cfg.BHhd == 1 && i < ien) {
                double sigW0 = (i < 7) ? sigW10d * pow((medtauds[i]+0.01)/10.01, betaW)
                                       : sigW0td;
                double hsigW = (i < 7) ? hsigWt : hsigWT;
                double sigvbs[3] = {};
                if (cfg.BHhb == 1 && i == 8) {
                    void calc_sigvb(double xb, double yb, double zb, double *sigvbs);
                    calc_sigvb(xyb[0], xyb[1], z, sigvbs);
                }
                double sigW     = (i == 8) ? sigvbs[2] : sigW0 * exp(-(R - R0) / hsigW);
                double sigW2    = sigW * sigW;
                double sigzkick2 = cfg.vkickBH * cfg.vkickBH * PI / 8;
                double sigzadd  = sqrt(sigW2 + sigzkick2);
                double fvBH     = pow(sigzadd / sigW, cfg.betaBH);
                double RhdBH    = (cfg.fixRhdBH == 1) ? cfg.RhdBH0
                                                       : cfg.RhdBH0 * (1 + sigW2 / sigzkick2);
                double fzdBH    = exp((R - R0) / RhdBH) * fvBH;
                double zd0 = (i < 8) ? zd[i] : 235.344943180979;
                double zdBH = (i < 8) ? zd0 * fzdBH : zd0 * fvBH;
                double rhoMS = 0, rhoBH = 0;
                if (i == 8) {
                    double x0E = 668.323640191308, y0E = 277.674592258175;
                    double C1E = 1.40903573470129, C2E = 3.3497118832179;
                    double xn = fabs(xyb[0] / x0E), yn = fabs(xyb[1] / y0E);
                    double Rs = pow(pow(xn, C1E) + pow(yn, C1E), 1.0/C1E);
                    double zn = fabs(z / zd0);
                    double rs = pow(pow(Rs, C2E) + pow(zn, C2E), 1.0/C2E);
                    rhoMS = exp(-rs);
                    zn = fabs(z / zdBH);
                    rs = pow(pow(Rs, C2E) + pow(zn, C2E), 1.0/C2E);
                    rhoBH = exp(-rs) / fvBH;
                } else {
                    rhoMS = (i < 7) ? 4.0 / (exp(2*z/zd0) + exp(-2*z/zd0) + 2)
                                    : exp(-fabs(z) / zd0);
                    rhoBH = (i < 7) ? 4.0 / (exp(2*z/zdBH) + exp(-2*z/zdBH) + 2)
                                    : exp(-fabs(z) / zdBH);
                    rhoBH /= fzdBH;
                }
                if (cfg.printBHfac)
                    printf("%d %5.0f %4.0f %5.2f %6.2f %6.2f %.6f %6.1f %6.1f"
                           " %.4e %.4e %.6f %.3f",
                           i, D_[ibin], R, z, sigW, sigzadd, fvBH, zd0, zdBH,
                           rhoMS, rhoBH, rhoBH/rhoMS, fBHtmp);
                fBHtmp *= rhoBH / rhoMS;
            }

            fBH_[i][ibin] = fBHtmp;
            if (i < 8 && cfg.printBHfac)
                printf(" %.6f\n", fBH_[i][ibin]);

            double nMS = (i == 8)  ? n0MSb  * rhos[8]
                       : (i == 9)  ? n0MSND * rhos[9]
                       : (i == 10) ? n0MSSH * rhos[10]
                       :             n0MSd[i] * rhos[i];
            double rho = (i == 8)  ? n0b   * rhos[8]
                       : (i == 9)  ? n0ND  * rhos[9]
                       : (i == 10) ? n0SH  * rhos[10]
                       :             n0d[i] * rhos[i];

            if (cfg.AI0 > 0 && cfg.Isen - cfg.Isst > 0 &&
                cfg.EVI0 > 0 && cfg.VIsen - cfg.VIsst > 0) {
                double fIVIs = fIVI_detect(extI, cfg.Isst, cfg.Isen,
                                           extVI, cfg.VIsst, cfg.VIsen, i);
                rhoD_S_[i][ibin] = nMS * fIVIs * 1e-06 * D_[ibin] * D_[ibin];
                nallS_ += rhoD_S_[i][ibin] * dD_;
            } else if (cfg.AI0 > 0 && cfg.Isen - cfg.Isst > 0) {
                double fIs = fLF_detect(extI, cfg.Isst, cfg.Isen, i);
                rhoD_S_[i][ibin] = nMS * fIs * 1e-06 * D_[ibin] * D_[ibin];
                nallS_ += rhoD_S_[i][ibin] * dD_;
            } else {
                double tmpDswt = (cfg.gammaDs == 0.5)
                    ? sqrt(D_[ibin] / 8000.0)
                    : pow((D_[ibin] + 10) / 8000.0, fabs(cfg.gammaDs));
                if (cfg.gammaDs < 0) tmpDswt = 1.0 / tmpDswt;
                rhoD_S_[i][ibin] = nMS * tmpDswt * 1e-03;
            }

            if (cfg.wtD_L != 0)
                rho *= pow((D_[ibin] + 1000) / 4500.0, cfg.wtD_L);
            rhoD_L_[i][ibin] = cfg.check_D
                ? rho * 1e-06 * D_[ibin] * D_[ibin]
                : rho;

            cumu_rho_S_[i][ibin] = (ibin == 0) ? 0.0
                : cumu_rho_S_[i][ibin-1] + 0.5*(rhoD_S_[i][ibin-1] + rhoD_S_[i][ibin])*dD_;
            cumu_rho_L_[i][ibin] = (ibin == 0) ? 0.0
                : cumu_rho_L_[i][ibin-1] + 0.5*(rhoD_L_[i][ibin-1] + rhoD_L_[i][ibin])*dD_;
            cumu_rho_all_S_[ibin] += cumu_rho_S_[i][ibin];
            cumu_rho_all_L_[ibin] += cumu_rho_L_[i][ibin];
            rhosum += cfg.printrhoS ? rhoD_S_[i][ibin] : rho;

            if (ibin % cfg.npri == 0) {
                if (cfg.printrhoS) {
                    printf(" %d: %.1e ", i, rhoD_S_[i][ibin]);
                    printf("( %.2e )", cumu_rho_S_[i][ibin]);
                } else {
                    printf(" %d: %.1e ", i, rho);
                    printf("( %.2e )", cumu_rho_L_[i][ibin]);
                }
            }
        }

        if (ibin % cfg.npri == 0) {
            printf(" All: %.1e ", rhosum);
            printf(cfg.printrhoS ? "( %.2e )\n" : "( %.2e )\n",
                   cfg.printrhoS ? cumu_rho_all_S_[ibin] : cumu_rho_all_L_[ibin]);
        }
    }

    if (cfg.printBHfac) exit(1);

    for (int i = 0; i < ncomp_; i++) {
        double norm_S = cumu_rho_S_[i][nbin_];
        double norm_L = cumu_rho_L_[i][nbin_];
        if (norm_S == 0 && i == 9)  continue;
        if (norm_S == 0 && i == 10) continue;
        for (int ibin = 0; ibin <= nbin_; ibin++) {
            double Pnorm_S = cumu_rho_S_[i][ibin] / norm_S;
            int intp_S = (int)(Pnorm_S * 20);
            if (ibinptiles_S_[i][intp_S] == 0)
                ibinptiles_S_[i][intp_S] = (intp_S == 0) ? 1 : (int)(ibin + 0.5);
            double Pnorm_L = cumu_rho_L_[i][ibin] / norm_L;
            int intp_L = (int)(Pnorm_L * 20);
            if (ibinptiles_L_[i][intp_L] == 0)
                ibinptiles_L_[i][intp_L] = (intp_L == 0) ? 1 : (int)(ibin + 0.5);
        }
    }
}

int LineOfSightDensityGrid::sample_source_component(double ran) const {
    double cumu = 0.0;
    int i_s;
    for (i_s = 0; i_s < ncomp_; i_s++) {
        cumu += cumu_rho_S_[i_s][nbin_] / cumu_rho_all_S_[nbin_];
        if (ran < cumu) break;
    }
    return i_s;
}

double LineOfSightDensityGrid::sample_source_distance(int i_s, double ran) const {
    double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
    int inttmp = (int)(ran * 20);
    int kst = 1;
    for (int itmp = inttmp; itmp > 0; itmp--) {
        kst = ibinptiles_S_[i_s][itmp];
        if (kst > 0) break;
    }
    double target = ran * cumu_rho_S_[i_s][nbin_];
    return getcumu2xist(nbin_ + 1,
                        const_cast<double*>(D_.data()),
                        const_cast<double*>(cumu_rho_S_[i_s].data()),
                        const_cast<double*>(rhoD_S_[i_s].data()),
                        target, kst, 0);
}

void LineOfSightDensityGrid::get_lens_importance_bins(
        double D_s, double thetaEmin, double thetaEmax,
        double piEmin, double piEmax,
        int &nbinDlmin, int &nbinDlmax) const {
    nbinDlmin = 0;
    nbinDlmax = (int)(D_s / dD_);
    if (thetaEmax - thetaEmin > 0 && piEmax - piEmin > 0) {
        double Dlmin = 1000.0 / (piEmax*thetaEmax + 1000.0/D_s);
        double Dlmax = 1000.0 / (piEmin*thetaEmin + 1000.0/D_s);
        nbinDlmin = (int)(Dlmin / dD_);
        nbinDlmax = (int)(Dlmax / dD_) + 1;
        int nbinDs = (int)(D_s / dD_);
        if (nbinDlmin < 0) nbinDlmin = 0;
        if (nbinDlmax > nbinDs) nbinDlmax = nbinDs;
    }
}

int LineOfSightDensityGrid::sample_lens_component(
        int nbinDlmin, int nbinDlmax, double ran) const {
    double total = cumu_rho_all_L_[nbinDlmax] - cumu_rho_all_L_[nbinDlmin];
    double cumu  = 0.0;
    int i_l;
    for (i_l = 0; i_l < ncomp_; i_l++) {
        cumu += (cumu_rho_L_[i_l][nbinDlmax] - cumu_rho_L_[i_l][nbinDlmin]) / total;
        if (ran < cumu) break;
    }
    return i_l;
}

double LineOfSightDensityGrid::sample_lens_distance(
        int i_l, int nbinDlmin, int nbinDlmax, double ran) const {
    double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
    double scaled = ran * (cumu_rho_L_[i_l][nbinDlmax] - cumu_rho_L_[i_l][nbinDlmin])
                   + cumu_rho_L_[i_l][nbinDlmin];
    int inttmp = (int)(scaled * 20.0 / cumu_rho_L_[i_l][nbin_]);
    int kst = 1;
    for (int itmp = inttmp; itmp > 0; itmp--) {
        kst = ibinptiles_L_[i_l][itmp];
        if (kst > 0) break;
    }
    return getcumu2xist(nbin_ + 1,
                        const_cast<double*>(D_.data()),
                        const_cast<double*>(cumu_rho_L_[i_l].data()),
                        const_cast<double*>(rhoD_L_[i_l].data()),
                        scaled, kst, 0);
}

double LineOfSightDensityGrid::importance_weight(int nbinDlmin, int nbinDlmax) const {
    return (cumu_rho_all_L_[nbinDlmax] - cumu_rho_all_L_[nbinDlmin])
         / cumu_rho_all_L_[nbin_];
}

double LineOfSightDensityGrid::fBH_at_distance(int comp, double D) const {
    double interp_x(int n, double *F, double xst, double dx, double xreq);
    return interp_x(nbin_ + 1,
                    const_cast<double*>(fBH_[comp].data()),
                    0.0, dD_, D);
}

} // namespace genulens
