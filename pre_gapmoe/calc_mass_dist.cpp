/* calc_mass_dist.cpp
 * Output the present-day mass function (PDMF) per Galactic component.
 * Integrates the IMF (broken power law) over the age distribution of each
 * component, applying the stellar turnoff mass from Minidie.dat.
 *
 * Usage (run from the genulens project root):
 *   calc_mass_dist [options]
 *
 * Key options:
 *   VERBOSITY <n>  0: PDMF table (default), 1: also print per-component summary
 *   Nmass <n>      Number of IMF log-mass intervals (default 1000)
 *   (plus all galaxy model options: model, DISK, addX, rhot0, alpha1..alpha4, ...)
 *
 * Output columns (tab-separated):
 *   logM[Msun]  dN/dlogM_d0..d6  dN/dlogM_thick  dN/dlogM_bar  dN/dlogM_NSD  dN/dlogM_halo  dN/dlogM_tot
 *
 * Units: [pc^-3] per unit dlogM, normalised to the solar position (R=R0, z=0).
 * Integrating each column over logM recovers n0MS for that component.
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
    printf("Output the present-day mass function per Galactic component.\n\n");
    printf("Key options:\n");
    printf("  VERBOSITY <n>  0: PDMF table, 1: add per-component summary (default 0)\n");
    printf("  Nmass <n>      Number of IMF log-mass intervals (default 1000)\n");
    printf("  -h, --help     Show this help and exit\n\n");
    printf("Useful galaxy model options include:\n");
    printf("  model <n>  DISK <n>  addX <n>  NSD <n>  SH <n>\n");
    printf("  M0/M1/M2/M3/Ml/Mu <Msun>  alpha1..alpha4 <slope>\n\n");
    printf("Resolution:\n");
    printf("  This tool is deterministic. It does not need AUTOERR; mass-grid accuracy is set by Nmass.\n\n");
    printf("Output columns:\n");
    printf("  logM[Msun]  dN/dlogM[0..10]  dN/dlogM_tot\n\n");
    printf("Component indices: 0-6 thin disk, 7 thick disk, 8 bar, 9 NSD, 10 halo.\n");
}

int main(int argc, char **argv)
{
    if (has_help_option(argc, argv)) {
        print_usage(argv[0]);
        return 0;
    }

    // ---- Initialize galaxy model (fills normalizations and Minidie tables) ----
    init_galactic_model(argc, argv, 0);

    int VERBOSE = (int)getOptiond(argc, argv, "VERBOSITY", 1, 0);

    // ---- Allocate PDMF array [ncomp][nm+1] ----
    double **pdmf = (double**)calloc(ncomp, sizeof(double*));
    for (int i = 0; i < ncomp; i++) pdmf[i] = (double*)calloc(nm+1, sizeof(double));

    // ---- Thin disk (components 0-6): exponential SFR, age bins ----
    {
        int iages7[7] = {15, 100, 200, 300, 500, 700, 1000};
        double gamma  = 1.0 / tSFR;
        for (int j = 1; j <= 1000; j++) {
            int itmp = (int)((j - agesD[0]) / (double)(agesD[1] - agesD[0]) + 0.5);
            if (itmp < 0)     itmp = 0;
            if (itmp >= nageD) itmp = nageD - 1;
            double logMdie = log10(MinidieD[itmp]);
            double wtSFR = exp(-gamma * (1000 - j) * 0.01);
            int idisk = (j<=iages7[0])?0:(j<=iages7[1])?1:(j<=iages7[2])?2
                       :(j<=iages7[3])?3:(j<=iages7[4])?4:(j<=iages7[5])?5:6;
            int iMdie = (int)((logMdie - logMst) / dlogM);
            if (iMdie > nm) iMdie = nm;
            if (iMdie < 0)  continue;
            for (int iM = 0; iM <= iMdie; iM++)
                pdmf[idisk][iM] += g_PlogM[iM] * wtSFR;
        }
        for (int i = 0; i < 7; i++) {
            double sum = 0;
            for (int iM = 0; iM <= nm; iM++) sum += pdmf[i][iM] * dlogM;
            double norm = (sum > 0) ? n0MSd[i] / sum : 0;
            for (int iM = 0; iM <= nm; iM++) pdmf[i][iM] *= norm;
        }
    }

    // ---- Thick disk (component 7): single old population ----
    {
        double logMdie = log10(MinidieD[nageD-2]);
        int iMdie = (int)((logMdie - logMst) / dlogM);
        if (iMdie > nm) iMdie = nm;
        for (int iM = 0; iM <= iMdie && iM <= nm; iM++) pdmf[7][iM] = g_PlogM[iM];
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[7][iM] * dlogM;
        double norm = (sum > 0) ? n0MSd[7] / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[7][iM] *= norm;
    }

    // ---- Bar (component 8): Gaussian SFR around mageB ----
    {
        for (int i = 0; i < nageB; i++) {
            double tau  = 0.01 * agesB[i];
            double dw   = (tau - mageB) / sageB;
            double wt   = exp(-0.5 * dw * dw);
            double logMdie = log10(MinidieB[i]);
            int iMdie = (int)((logMdie - logMst) / dlogM);
            if (iMdie > nm) iMdie = nm;
            if (iMdie < 0)  continue;
            for (int iM = 0; iM <= iMdie; iM++)
                pdmf[8][iM] += g_PlogM[iM] * wt;
        }
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[8][iM] * dlogM;
        double norm = (sum > 0) ? n0MSb / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[8][iM] *= norm;
    }

    // ---- NSD (component 9): Gaussian SFR around mageND ----
    if (ND > 0 && n0MSND > 0) {
        for (int i = 0; i < nageND; i++) {
            double tau  = 0.01 * agesND[i];
            double dw   = (tau - mageND) / sageND;
            double wt   = exp(-0.5 * dw * dw);
            double logMdie = log10(MinidieND[i]);
            int iMdie = (int)((logMdie - logMst) / dlogM);
            if (iMdie > nm) iMdie = nm;
            if (iMdie < 0)  continue;
            for (int iM = 0; iM <= iMdie; iM++)
                pdmf[9][iM] += g_PlogM[iM] * wt;
        }
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[9][iM] * dlogM;
        double norm = (sum > 0) ? n0MSND / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[9][iM] *= norm;
    }

    // ---- Stellar halo (component 10): single very old population ----
    if (SH > 0) {
        double logMdie = log10(MinidieD[nageD-1]);
        int iMdie = (int)((logMdie - logMst) / dlogM);
        if (iMdie > nm) iMdie = nm;
        for (int iM = 0; iM <= iMdie && iM <= nm; iM++) pdmf[10][iM] = g_PlogM[iM];
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[10][iM] * dlogM;
        double norm = (sum > 0) ? n0MSSH / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[10][iM] *= norm;
    }

    // ---- Header ----
    printf("# calc_mass_dist\n");
    printf("# DISK=%d  hDISK=%d  model=%d  addX=%d  ND=%d  SH=%d\n",
           DISK, hDISK, model, addX, ND, SH);
    printf("# IMF: logM in [%.4f, %.4f] with dlogM=%.5f (%d bins)\n",
           logMst, logMst + nm*dlogM, dlogM, nm);

    if (VERBOSE >= 1) {
        // Per-component summary: n0MS, aveMMS, n0WD, n0RG
        const char *names[11] = {
            "disk0","disk1","disk2","disk3","disk4","disk5","disk6",
            "thick","bar","NSD","halo"
        };
        double n0WD[11], n0RG[11], n0MS_check[11], aveMMS[11];
        // n0WD from globals
        for (int i = 0; i < 7; i++)  { n0WD[i]=n0d[i]-n0MSd[i]; n0RG[i]=n0RGd[i]; }
        n0WD[7]=n0d[7]-n0MSd[7]; n0RG[7]=n0RGd[7];
        n0WD[8]=n0b -n0MSb;  n0RG[8]=n0RGb;
        n0WD[9]=n0ND-n0MSND; n0RG[9]=n0RGND;
        n0WD[10]=n0SH-n0MSSH; n0RG[10]=n0RGSH;
        // aveMMS from PDMF
        double n0MS_arr[11] = {n0MSd[0],n0MSd[1],n0MSd[2],n0MSd[3],n0MSd[4],n0MSd[5],n0MSd[6],
                               n0MSd[7],n0MSb,n0MSND,n0MSSH};
        for (int i = 0; i < ncomp; i++) {
            double sumN = 0, sumMN = 0;
            for (int iM = 0; iM <= nm; iM++) {
                double M = pow(10.0, g_logMass[iM]);
                sumN  += pdmf[i][iM] * dlogM;
                sumMN += M * pdmf[i][iM] * dlogM;
            }
            n0MS_check[i] = sumN;
            aveMMS[i] = (sumN > 0) ? sumMN / sumN : 0;
        }
        printf("# %-6s  %12s  %12s  %12s  %12s  %12s\n",
               "comp","n0MS[pc-3]","aveMMS[Msun]","n0WD[pc-3]","n0RG[pc-3]","n0tot[pc-3]");
        for (int i = 0; i < ncomp; i++) {
            double n0tot = n0MS_arr[i] + n0WD[i];
            printf("# %-6s  %12.4e  %12.4f  %12.4e  %12.4e  %12.4e\n",
                   names[i], n0MS_arr[i], aveMMS[i], n0WD[i], n0RG[i], n0tot);
        }
    }

    // ---- Column header ----
    printf("# Columns:\n# logM[Msun]");
    for (int i = 0; i < ncomp; i++) printf("  dN/dlogM[%d]", i);
    printf("  dN/dlogM_tot\n");

    // ---- Main output ----
    for (int iM = 0; iM <= nm; iM++) {
        double tot = 0;
        printf("%.6f", g_logMass[iM]);
        for (int i = 0; i < ncomp; i++) {
            printf("\t%.6e", pdmf[i][iM]);
            tot += pdmf[i][iM];
        }
        printf("\t%.6e\n", tot);
    }

    // ---- Cleanup ----
    for (int i = 0; i < ncomp; i++) free(pdmf[i]);
    free(pdmf);
    free(g_logMass); free(g_PlogM); free(g_PlogM_cum_norm);
    gsl_rng_free(r_rng);
    return 0;
}
