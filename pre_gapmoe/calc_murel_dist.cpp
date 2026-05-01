/* calc_murel_dist.cpp
 * Output the relative proper-motion (murel) distribution over a Dl × Ds grid
 * (default) or at a single (Dl, Ds) pair, sampled via Monte Carlo using the
 * Shu distribution function kinematics.
 *
 * Usage (run from the genulens project root):
 *   calc_murel_dist [options]
 *
 * Key options:
 *   l <deg>          Galactic longitude (default 1.0)
 *   b <deg>          Galactic latitude  (default -3.9)
 *   Nsimu <n>        MC draws per (Dl, Ds) cell (default 10000000)
 *   mumax <mas/yr>   Upper edge of murel histogram (default 300)
 *   dmu  <mas/yr>    Bin width (default 0.5)
 *   GRID <n>         1: Dl×Ds grid mode (default), 0: single (Dl, Ds) point
 *
 * Grid mode options (GRID=1):
 *   DLmin <pc>       Lens dist grid start    (default 0)
 *   DLmax <pc>       Lens dist grid end      (default 12000)
 *   DLstep <pc>      Lens dist bin size      (default 500)
 *   DSmin <pc>       Source dist grid start  (default 0)
 *   DSmax <pc>       Source dist grid end    (default 16000)
 *   DSstep <pc>      Source dist bin size    (default 500)
 *
 * Single-point options (GRID=0):
 *   Ds <pc>          Source distance (default 8000)
 *   Dl <pc>          Lens   distance (default 4000)
 *   VERBOSITY <n>    0: murel histogram (default)
 *                    1: + lens component weights + per-component breakdown
 *                    2: raw event list (wt murel mul mub phi i_l i_s tau_l tau_s)
 *
 * (plus all galaxy model options: model, DISK, addX, ...)
 *
 * Output (GRID=1):
 *   DS[pc]  DL[pc]  murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi
 *   If AUTOERR=1, append: mu_count  mu_relerr  phi_count  phi_relerr  Ndraw
 *
 * Output (GRID=0, VERBOSITY=0 or 1):
 *   murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi  [per-component columns if VERBOSITY=1]
 *   If AUTOERR=1, append: mu_count  mu_relerr  phi_count  phi_relerr  Ndraw
 *
 * Output (GRID=0, VERBOSITY=2):
 *   wt  murel  mu_l  mu_b  phi  i_l  i_s  tau_l  tau_s
 */

#include "galactic_model.h"
#include "option.h"

static int has_help_option(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "help"))
            return 1;
    }
    return 0;
}

static void print_usage(const char *prog)
{
    printf("Usage: %s [options]\n\n", prog);
    printf("Output the relative proper-motion distribution over a Dl x Ds grid or a single point.\n\n");
    printf("Key options:\n");
    printf("  l <deg>          Galactic longitude (default 1.0)\n");
    printf("  b <deg>          Galactic latitude  (default -3.9)\n");
    printf("  Nsimu <n>        MC draws per cell or point (default 10000000)\n");
    printf("  mumax <mas/yr>   Upper edge of murel histogram (default 300)\n");
    printf("  dmu <mas/yr>     Bin width (default 0.5)\n");
    printf("  GRID <n>         1: Dl x Ds grid, 0: single point (default 1)\n");
    printf("  -h, --help       Show this help and exit\n\n");
    printf("Automatic precision options:\n");
    printf("  AUTOERR <n>      1: stop early once histogram precision is reached (default 1)\n");
    printf("  ERR_TARGET <x>   Target relative Poisson error per active bin (default 0.03)\n");
    printf("  ERR_CHECK <n>    Check precision every n accepted draws (default 100000)\n");
    printf("  ERR_CDFLO <x>    Ignore lower CDF tail for precision checks (default 0.01)\n");
    printf("  ERR_CDFHI <x>    Ignore upper CDF tail for precision checks (default 0.99)\n");
    printf("  ERR_MINPROB <x>  Ignore bins with probability below this value (default 1e-4)\n");
    printf("  ERR_FRAC <x>     Required fraction of active bins passing target (default 0.90)\n\n");
    printf("Grid mode options (GRID=1):\n");
    printf("  DLmin/DLmax/DLstep <pc>  Lens grid, defaults 0/12000/500\n");
    printf("  DSmin/DSmax/DSstep <pc>  Source grid, defaults 0/16000/500\n\n");
    printf("Single-point options (GRID=0):\n");
    printf("  Dl <pc>          Lens distance (default 4000)\n");
    printf("  Ds <pc>          Source distance (default 8000)\n");
    printf("  VERBOSITY <n>    0: histogram, 1: component columns, 2: raw events\n\n");
    printf("Output columns:\n");
    printf("  GRID=1: DS[pc]  DL[pc]  murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi\n");
    printf("  GRID=0: murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi  [component columns]\n");
    printf("  AUTOERR=1 appends: mu_count  mu_relerr  phi_count  phi_relerr  Ndraw\n");
    printf("  VERBOSITY=2: wt  murel  mu_l  mu_b  phi  i_l  i_s  tau_l  tau_s\n\n");
    printf("AUTOERR requires both the murel and phi histograms to pass the same precision target.\n");
    printf("Precision comments include mu_overflow for draws outside the murel histogram range.\n");
    printf("murel and phi are separate 1D histograms printed side by side; rows beyond one histogram use 0 for that side.\n");
    printf("phi is atan2(mu_E, mu_N), using the same PA rotation as murel_sampling.c.\n");
}

static int check_hist_precision(double *hist, int nbins, double total_draws,
                                double target_relerr, double cdf_lo, double cdf_hi,
                                double min_prob, double required_frac,
                                double *max_relerr_out, double *pass_frac_out,
                                int *n_active_out)
{
    double hist_sum = 0.0;
    for (int ib = 0; ib < nbins; ib++) hist_sum += hist[ib];
    if (hist_sum <= 0 || total_draws <= 0) {
        *max_relerr_out = -1.0;
        *pass_frac_out = 0.0;
        *n_active_out = 0;
        return 0;
    }

    int n_active = 0, n_pass = 0;
    double max_relerr = 0.0;
    double cdf = 0.0;
    for (int ib = 0; ib < nbins; ib++) {
        double cdf_mid = (cdf + 0.5 * hist[ib]) / hist_sum;
        cdf += hist[ib];
        if (cdf_mid < cdf_lo || cdf_mid > cdf_hi) continue;

        double p = hist[ib] / total_draws;
        if (p < min_prob || hist[ib] <= 0) continue;

        double relerr = 1.0 / sqrt(hist[ib]);
        if (relerr > max_relerr) max_relerr = relerr;
        n_active++;
        if (relerr <= target_relerr) n_pass++;
    }

    *max_relerr_out = (n_active > 0) ? max_relerr : -1.0;
    *pass_frac_out = (n_active > 0) ? (double)n_pass / n_active : 0.0;
    *n_active_out = n_active;
    return (n_active > 0 && *pass_frac_out >= required_frac);
}

int main(int argc, char **argv)
{
    if (has_help_option(argc, argv)) {
        print_usage(argv[0]);
        return 0;
    }

    init_galactic_model(argc, argv, 1);

    double lSIMU  = getOptiond(argc,argv,"l",         1, 1.0);
    double bSIMU  = getOptiond(argc,argv,"b",         1, -3.9);
    long   Nsimu  = (long)getOptiond(argc,argv,"Nsimu",    1, 10000000);
    int    VERBOSE= (int)getOptiond(argc,argv,"VERBOSITY", 1, 0);
    double mumax  = getOptiond(argc,argv,"mumax",     1, 300.0);
    double dmu    = getOptiond(argc,argv,"dmu",       1, 0.5);
    int    GRID   = (int)getOptiond(argc,argv,"GRID",  1, 1);
    int    AUTOERR = getOptioni(argc,argv,"AUTOERR", 1, 1);
    double ERR_TARGET  = getOptiond(argc,argv,"ERR_TARGET",  1, 0.03);
    long   ERR_CHECK   = (long)getOptiond(argc,argv,"ERR_CHECK", 1, 100000);
    double ERR_CDFLO   = getOptiond(argc,argv,"ERR_CDFLO",   1, 0.01);
    double ERR_CDFHI   = getOptiond(argc,argv,"ERR_CDFHI",   1, 0.99);
    double ERR_MINPROB = getOptiond(argc,argv,"ERR_MINPROB", 1, 1e-4);
    double ERR_FRAC    = getOptiond(argc,argv,"ERR_FRAC",    1, 0.90);
    if (ERR_CHECK <= 0) ERR_CHECK = 100000;

    // Grid mode parameters
    double DLmin  = getOptiond(argc,argv,"DLmin",  1, 0.0);
    double DLmax  = getOptiond(argc,argv,"DLmax",  1, 12000.0);
    double DLstep = getOptiond(argc,argv,"DLstep", 1, 500.0);
    double DSmin  = getOptiond(argc,argv,"DSmin",  1, 0.0);
    double DSmax  = getOptiond(argc,argv,"DSmax",  1, 16000.0);
    double DSstep = getOptiond(argc,argv,"DSstep", 1, 500.0);

    // Single-point parameters
    double Ds     = getOptiond(argc,argv,"Ds", 1, 8000.0);
    double Dl     = getOptiond(argc,argv,"Dl", 1, 4000.0);

    int nbins = (int)(mumax / dmu) + 1;

    // Age lookup: median age per component (Gyr)
    double tau_median[11];
    for (int i = 0; i < 8; i++) tau_median[i] = medtauds[i];
    tau_median[8]  = mageB;
    tau_median[9]  = mageND;
    tau_median[10] = 12.0;

    double n0MS_arr[11] = {
        n0MSd[0],n0MSd[1],n0MSd[2],n0MSd[3],n0MSd[4],n0MSd[5],n0MSd[6],
        n0MSd[7], n0MSb, n0MSND, n0MSSH
    };

    double cosl = cos(lSIMU/180.0*PI), sinl = sin(lSIMU/180.0*PI);
    double cosb = cos(bSIMU/180.0*PI), sinb = sin(bSIMU/180.0*PI);
    double PA = 0.0, cosPA = 1.0, sinPA = 0.0;
    calc_PA(lSIMU, bSIMU, &PA, &cosPA, &sinPA);
    double phimin = -PI, phimax = PI, dphi = 5.0 * PI / 180.0;
    int nphibins = (int)((phimax - phimin) / dphi);

    // ---- Common header ----
    printf("# calc_murel_dist\n");
    printf("# (l, b) = (%.4f, %.4f) deg\n", lSIMU, bSIMU);
    printf("# Nsimu = %ld per cell, dmu = %.2f mas/yr\n", Nsimu, dmu);
    printf("# phi bins: [%.6f, %.6f] rad, dphi = %.6f rad (%d bins), PA = %.6f deg\n",
           phimin, phimax, dphi, nphibins, PA);
    printf("# DISK=%d  hDISK=%d  model=%d  addX=%d  ND=%d  SH=%d\n",
           DISK, hDISK, model, addX, ND, SH);
    if (AUTOERR) {
        printf("# AUTOERR=1: target_relerr=%.3f, check_every=%ld, CDF=[%.3f, %.3f], min_prob=%.2e, required_frac=%.3f\n",
               ERR_TARGET, ERR_CHECK, ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC);
        printf("# Nsimu is treated as the maximum accepted draw count.\n");
    }

    // ====================================================================
    // GRID MODE (default)
    // ====================================================================
    if (GRID == 1) {
        int nDL = (int)((DLmax - DLmin) / DLstep);
        int nDS = (int)((DSmax - DSmin) / DSstep);

        printf("# Grid: DL [%.0f, %.0f] step %.0f pc (%d bins)\n",
               DLmin, DLmax, DLstep, nDL);
        printf("# Grid: DS [%.0f, %.0f] step %.0f pc (%d bins)\n",
               DSmin, DSmax, DSstep, nDS);
        printf("# Columns:\n# DS[pc]  DL[pc]  murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi");
        if (AUTOERR) printf("  mu_count  mu_relerr  phi_count  phi_relerr  Ndraw");
        printf("\n");

        double *hist = (double*)calloc(nbins, sizeof(double));
        double *phi_hist = (double*)calloc(nphibins, sizeof(double));

        for (int j_DS = 0; j_DS < nDS; j_DS++) {
            double DS_c = DSmin + (j_DS + 0.5) * DSstep;
        for (int i_DL = 0; i_DL < nDL; i_DL++) {
            double DL_c = DLmin + (i_DL + 0.5) * DLstep;
            if (DL_c >= DS_c) continue;

            // Component weights at bin centers
            double *rhos_L = (double*)calloc(ncomp, sizeof(double));
            double *rhos_S = (double*)calloc(ncomp, sizeof(double));
            double xyz[3]={}, xyb[2]={};
            calc_rho_each(DL_c, 0, rhos_L, xyz, xyb);
            calc_rho_each(DS_c, 0, rhos_S, xyz, xyb);

            double wL[11]={}, wS[11]={}, wL_tot=0, wS_tot=0;
            for (int i = 0; i < ncomp; i++) {
                wL[i] = n0MS_arr[i] * rhos_L[i];
                wS[i] = n0MS_arr[i] * rhos_S[i];
                wL_tot += wL[i];
                wS_tot += wS[i];
            }
            free(rhos_L); free(rhos_S);
            if (wL_tot <= 0 || wS_tot <= 0) continue;

            double cumwL[11]={}, cumwS[11]={};
            for (int i = 0; i < ncomp; i++) {
                cumwL[i] = (i==0 ? 0 : cumwL[i-1]) + wL[i]/wL_tot;
                cumwS[i] = (i==0 ? 0 : cumwS[i-1]) + wS[i]/wS_tot;
            }

            for (int ib = 0; ib < nbins; ib++) hist[ib] = 0;
            for (int ib = 0; ib < nphibins; ib++) phi_hist[ib] = 0;
            double wtot = 0;
            double mu_overflow = 0.0;

            double mu_max_relerr = 1e99, mu_pass_frac = 0.0;
            double phi_max_relerr = 1e99, phi_pass_frac = 0.0;
            int mu_n_active = 0, phi_n_active = 0;
            int converged_mu = 0, converged_phi = 0, converged = 0;

            while (wtot < Nsimu) {
                double ran_L = ran1();
                int i_l = ncomp - 1;
                for (int i = 0; i < ncomp; i++) if (ran_L < cumwL[i]) { i_l = i; break; }

                double ran_S = ran1();
                int i_s = ncomp - 1;
                for (int i = 0; i < ncomp; i++) if (ran_S < cumwS[i]) { i_s = i; break; }

                if ((i_l == 9 || i_s == 9) && ND == 0) continue;

                double vxyz_L[3]={}, vxyz_S[3]={};
                get_vxyz_ran(vxyz_L, i_l, tau_median[i_l], DL_c, lSIMU, bSIMU);
                get_vxyz_ran(vxyz_S, i_s, tau_median[i_s], DS_c, lSIMU, bSIMU);

                double dvx = vxyz_L[0] - vxyz_S[0];
                double dvy = vxyz_L[1] - vxyz_S[1];
                double dvz = vxyz_L[2] - vxyz_S[2];

                double vel_l = (-sinl)*dvx + cosl*dvy;
                double vel_b = (-sinb*cosl)*dvx + (-sinb*sinl)*dvy + cosb*dvz;

                double mul = vel_l * KS2MY / DL_c;
                double mub = vel_b * KS2MY / DL_c;
                double murel = sqrt(mul*mul + mub*mub);
                double muN = mub * cosPA + mul * sinPA;
                double muE = -mub * sinPA + mul * cosPA;
                double phi = atan2(muE, muN);

                wtot += 1.0;
                int ibin = (int)(murel / dmu);
                if (ibin < nbins) hist[ibin] += 1.0;
                else mu_overflow += 1.0;
                int iphi = (int)((phi - phimin) / dphi);
                if (iphi >= 0 && iphi < nphibins) phi_hist[iphi] += 1.0;

                if (AUTOERR && ((long)wtot % ERR_CHECK == 0)) {
                    converged_mu = check_hist_precision(hist, nbins, wtot, ERR_TARGET,
                                                        ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                        &mu_max_relerr, &mu_pass_frac, &mu_n_active);
                    converged_phi = check_hist_precision(phi_hist, nphibins, wtot, ERR_TARGET,
                                                         ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                         &phi_max_relerr, &phi_pass_frac, &phi_n_active);
                    converged = converged_mu && converged_phi;
                    if (converged) break;
                }
            }

            if (wtot <= 0) continue;
            if (AUTOERR && !converged) {
                converged_mu = check_hist_precision(hist, nbins, wtot, ERR_TARGET,
                                                    ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                    &mu_max_relerr, &mu_pass_frac, &mu_n_active);
                converged_phi = check_hist_precision(phi_hist, nphibins, wtot, ERR_TARGET,
                                                     ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                     &phi_max_relerr, &phi_pass_frac, &phi_n_active);
                converged = converged_mu && converged_phi;
            }
            if (AUTOERR) {
                printf("# precision DS=%.1f DL=%.1f Ndraw=%.0f "
                       "mu_overflow=%.0f mu_overflow_frac=%.6f "
                       "mu_active_bins=%d mu_max_relerr=%.4f mu_pass_frac=%.3f mu_converged=%d "
                       "phi_active_bins=%d phi_max_relerr=%.4f phi_pass_frac=%.3f phi_converged=%d "
                       "converged=%d\n",
                       DS_c, DL_c, wtot,
                       mu_overflow, mu_overflow / wtot,
                       mu_n_active, mu_max_relerr, mu_pass_frac, converged_mu,
                       phi_n_active, phi_max_relerr, phi_pass_frac, converged_phi,
                       converged);
            }
            for (int ib = 0; ib < nbins; ib++) hist[ib] /= (wtot * dmu);
            for (int ib = 0; ib < nphibins; ib++) phi_hist[ib] /= (wtot * dphi);

            int nout = (nbins > nphibins) ? nbins : nphibins;
            for (int ib = 0; ib < nout; ib++) {
                double mu_center = (ib < nbins) ? (ib + 0.5) * dmu : 0.0;
                double mu_density = (ib < nbins) ? hist[ib] : 0.0;
                double phi_center = (ib < nphibins) ? phimin + (ib + 0.5) * dphi : 0.0;
                double phi_density = (ib < nphibins) ? phi_hist[ib] : 0.0;
                printf("%.1f\t%.1f\t%.3f\t%.6f\t%.6e\t%.6e",
                       DS_c, DL_c, mu_center, phi_center, mu_density, phi_density);
                if (AUTOERR) {
                    double mu_count = (ib < nbins) ? hist[ib] * wtot * dmu : 0.0;
                    double mu_relerr = (mu_count > 0) ? 1.0 / sqrt(mu_count) : 0.0;
                    double phi_count = (ib < nphibins) ? phi_hist[ib] * wtot * dphi : 0.0;
                    double phi_relerr = (phi_count > 0) ? 1.0 / sqrt(phi_count) : 0.0;
                    printf("\t%.0f\t%.6e\t%.0f\t%.6e\t%.0f",
                           mu_count, mu_relerr, phi_count, phi_relerr, wtot);
                }
                printf("\n");
            }
        }} // i_DL, j_DS

        free(hist);
        free(phi_hist);

    // ====================================================================
    // SINGLE-POINT MODE (GRID=0)
    // ====================================================================
    } else {
        if (Dl >= Ds) { printf("# ERROR: Dl must be < Ds\n"); return 1; }

        double pirel = (1.0/Dl - 1.0/Ds) * 1000.0;
        printf("# Ds = %.0f pc, Dl = %.0f pc, pirel = %.4f mas\n", Ds, Dl, pirel);

        double *rhos_L = (double*)calloc(ncomp, sizeof(double));
        double *rhos_S = (double*)calloc(ncomp, sizeof(double));
        double xyz_L[3]={}, xyz_S[3]={}, xyb_L[2]={}, xyb_S[2]={};
        calc_rho_each(Dl, 0, rhos_L, xyz_L, xyb_L);
        calc_rho_each(Ds, 0, rhos_S, xyz_S, xyb_S);

        double wL[11]={}, wS[11]={}, wL_tot=0, wS_tot=0;
        for (int i = 0; i < ncomp; i++) {
            wL[i] = n0MS_arr[i] * rhos_L[i];
            wS[i] = n0MS_arr[i] * rhos_S[i];
            wL_tot += wL[i];
            wS_tot += wS[i];
        }
        free(rhos_L); free(rhos_S);
        if (wL_tot <= 0 || wS_tot <= 0) {
            printf("# ERROR: zero density at Dl=%.0f or Ds=%.0f pc\n", Dl, Ds);
            return 1;
        }

        double cumwL[11]={}, cumwS[11]={};
        for (int i = 0; i < ncomp; i++) {
            cumwL[i] = (i==0 ? 0 : cumwL[i-1]) + wL[i]/wL_tot;
            cumwS[i] = (i==0 ? 0 : cumwS[i-1]) + wS[i]/wS_tot;
        }

        double *hist    = (double*)calloc(nbins, sizeof(double));
        double *hist_wt = (double*)calloc(nbins, sizeof(double));
        double *phi_hist = (double*)calloc(nphibins, sizeof(double));
        double **hist_comp = NULL;
        if (VERBOSE >= 1) {
            hist_comp = (double**)calloc(ncomp, sizeof(double*));
            for (int i = 0; i < ncomp; i++)
                hist_comp[i] = (double*)calloc(nbins, sizeof(double));
        }
        double wtot = 0;
        double mu_overflow = 0.0;
        double mu_max_relerr = 1e99, mu_pass_frac = 0.0;
        double phi_max_relerr = 1e99, phi_pass_frac = 0.0;
        int mu_n_active = 0, phi_n_active = 0;
        int converged_mu = 0, converged_phi = 0, converged = 0;

        while (wtot < Nsimu) {
            double ran_L = ran1();
            int i_l = ncomp - 1;
            for (int i = 0; i < ncomp; i++) if (ran_L < cumwL[i]) { i_l = i; break; }

            double ran_S = ran1();
            int i_s = ncomp - 1;
            for (int i = 0; i < ncomp; i++) if (ran_S < cumwS[i]) { i_s = i; break; }

            if ((i_l == 9 || i_s == 9) && ND == 0) continue;

            double tau_l = tau_median[i_l];
            double tau_s = tau_median[i_s];

            double vxyz_L[3]={}, vxyz_S[3]={};
            get_vxyz_ran(vxyz_L, i_l, tau_l, Dl, lSIMU, bSIMU);
            get_vxyz_ran(vxyz_S, i_s, tau_s, Ds, lSIMU, bSIMU);

            double dvx = vxyz_L[0] - vxyz_S[0];
            double dvy = vxyz_L[1] - vxyz_S[1];
            double dvz = vxyz_L[2] - vxyz_S[2];

            double vel_l = (-sinl)*dvx + cosl*dvy;
            double vel_b = (-sinb*cosl)*dvx + (-sinb*sinl)*dvy + cosb*dvz;

            double mul = vel_l * KS2MY / Dl;
            double mub = vel_b * KS2MY / Dl;
            double murel = sqrt(mul*mul + mub*mub);
            double muN = mub * cosPA + mul * sinPA;
            double muE = -mub * sinPA + mul * cosPA;
            double phi = atan2(muE, muN);

            double wt = 1.0;
            wtot += wt;

            if (VERBOSE >= 2) {
                printf("%.6e %.4f %.4f %.4f %.6f %d %d %.3f %.3f\n",
                       wt, murel, mul, mub, phi, i_l, i_s, tau_l, tau_s);
                continue;
            }

            int ibin = (int)(murel / dmu);
            if (ibin < nbins) {
                hist[ibin]    += wt;
                hist_wt[ibin] += wt;
                if (VERBOSE >= 1) hist_comp[i_l][ibin] += wt;
            } else mu_overflow += wt;
            int iphi = (int)((phi - phimin) / dphi);
            if (iphi >= 0 && iphi < nphibins) phi_hist[iphi] += wt;

            if (AUTOERR && VERBOSE < 2 && ((long)wtot % ERR_CHECK == 0)) {
                converged_mu = check_hist_precision(hist, nbins, wtot, ERR_TARGET,
                                                    ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                    &mu_max_relerr, &mu_pass_frac, &mu_n_active);
                converged_phi = check_hist_precision(phi_hist, nphibins, wtot, ERR_TARGET,
                                                     ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                     &phi_max_relerr, &phi_pass_frac, &phi_n_active);
                converged = converged_mu && converged_phi;
                if (converged) break;
            }
        }

        if (VERBOSE >= 2) {
            free(hist); free(hist_wt); free(phi_hist);
            free(g_logMass); free(g_PlogM); free(g_PlogM_cum_norm);
            gsl_rng_free(r_rng);
            return 0;
        }

        if (AUTOERR) {
            if (!converged) {
                converged_mu = check_hist_precision(hist, nbins, wtot, ERR_TARGET,
                                                    ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                    &mu_max_relerr, &mu_pass_frac, &mu_n_active);
                converged_phi = check_hist_precision(phi_hist, nphibins, wtot, ERR_TARGET,
                                                     ERR_CDFLO, ERR_CDFHI, ERR_MINPROB, ERR_FRAC,
                                                     &phi_max_relerr, &phi_pass_frac, &phi_n_active);
                converged = converged_mu && converged_phi;
            }
            printf("# precision Ndraw=%.0f "
                   "mu_overflow=%.0f mu_overflow_frac=%.6f "
                   "mu_active_bins=%d mu_max_relerr=%.4f mu_pass_frac=%.3f mu_converged=%d "
                   "phi_active_bins=%d phi_max_relerr=%.4f phi_pass_frac=%.3f phi_converged=%d "
                   "converged=%d\n",
                   wtot,
                   mu_overflow, mu_overflow / wtot,
                   mu_n_active, mu_max_relerr, mu_pass_frac, converged_mu,
                   phi_n_active, phi_max_relerr, phi_pass_frac, converged_phi,
                   converged);
        }

        for (int ib = 0; ib < nbins; ib++) {
            hist[ib] /= (wtot * dmu);
            if (VERBOSE >= 1)
                for (int i = 0; i < ncomp; i++) hist_comp[i][ib] /= (wtot * dmu);
        }
        for (int ib = 0; ib < nphibins; ib++) phi_hist[ib] /= (wtot * dphi);

        if (VERBOSE >= 1) {
            printf("# Lens component weights at Dl:\n");
            for (int i = 0; i < ncomp; i++)
                printf("#   comp%2d: %.4e\n", i, wL[i]/wL_tot);
        }

        printf("# Columns:\n# murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi");
        if (VERBOSE >= 1)
            for (int i = 0; i < ncomp; i++) printf("  dP/dmurel[L=%d]", i);
        if (AUTOERR) printf("  mu_count  mu_relerr  phi_count  phi_relerr  Ndraw");
        printf("\n");

        int nout = (nbins > nphibins) ? nbins : nphibins;
        for (int ib = 0; ib < nout; ib++) {
            double mu_center = (ib < nbins) ? (ib + 0.5) * dmu : 0.0;
            double mu_density = (ib < nbins) ? hist[ib] : 0.0;
            double phi_center = (ib < nphibins) ? phimin + (ib + 0.5) * dphi : 0.0;
            double phi_density = (ib < nphibins) ? phi_hist[ib] : 0.0;
            printf("%.3f\t%.6f\t%.6e\t%.6e", mu_center, phi_center, mu_density, phi_density);
            if (VERBOSE >= 1)
                for (int i = 0; i < ncomp; i++) {
                    double comp_density = (ib < nbins) ? hist_comp[i][ib] : 0.0;
                    printf("\t%.6e", comp_density);
                }
            if (AUTOERR) {
                double mu_count = (ib < nbins) ? hist[ib] * wtot * dmu : 0.0;
                double mu_relerr = (mu_count > 0) ? 1.0 / sqrt(mu_count) : 0.0;
                double phi_count = (ib < nphibins) ? phi_hist[ib] * wtot * dphi : 0.0;
                double phi_relerr = (phi_count > 0) ? 1.0 / sqrt(phi_count) : 0.0;
                printf("\t%.0f\t%.6e\t%.0f\t%.6e\t%.0f",
                       mu_count, mu_relerr, phi_count, phi_relerr, wtot);
            }
            printf("\n");
        }

        free(hist); free(hist_wt);
        free(phi_hist);
        if (VERBOSE >= 1) {
            for (int i = 0; i < ncomp; i++) free(hist_comp[i]);
            free(hist_comp);
        }
    }

    free(g_logMass); free(g_PlogM); free(g_PlogM_cum_norm);
    gsl_rng_free(r_rng);
    return 0;
}
