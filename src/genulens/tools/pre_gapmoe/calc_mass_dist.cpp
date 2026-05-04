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

#include "genulens/cli/option.h"
#include "genulens/io/input_data.hpp"
#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/internal/runtime.hpp"

#define fopen(path, mode) genulens::open_input_file((path), (mode))

using namespace genulens;

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

    // ---- Initialize galaxy model ----
    RunContext ctx = Initializer().create_context();
    Initializer().initialize_rng(ctx, argc, argv);
    Initializer().read_model_options(ctx, argc, argv);
    Initializer().finalize_spatial_model(ctx, argc, argv);

    // lDs/bDs for density calculation (single sightline, not used by PDMF but required by init)
    ctx.density.lDs = (double*)malloc(sizeof(double));
    ctx.density.bDs = (double*)malloc(sizeof(double));
    ctx.density.lDs[0] = getOptiond(argc, argv, "l", 1, 1.0);
    ctx.density.bDs[0] = getOptiond(argc, argv, "b", 1, -3.9);

    // NSD option
    ctx.density.ND = getOptiond(argc, argv, "NSD", 1, 0);
    if (ctx.density.ND > 0) ctx.density.ND = 3;
    double MND = (ctx.density.ND == 1) ? 2.0e9 : (ctx.density.ND == 2) ? 7.0e8 : 0.0;
    MND = getOptiond(argc, argv, "MND", 1, MND);
    if (ctx.density.ND)
        ctx.density.rho0ND = (ctx.density.ND == 3) ? 1.0 : 0.25*MND/3.1415926535897932385/ctx.density.x0ND/ctx.density.y0ND/ctx.density.z0ND;

    PopulationRuntime population;
    population.initialize_mass_function(ctx, ctx.imf_options);
    Initializer().finalize_density_normalization(ctx);

    int VERBOSE = (int)getOptiond(argc, argv, "VERBOSITY", 1, 0);
    int ncomp = ctx.density.ncomp;
    int nm    = ctx.stellar.nm;
    double logMst = ctx.stellar.logMst;
    double dlogM  = ctx.stellar.dlogM;
    double *g_logMass        = population.log_mass;
    double *g_PlogM          = population.mass_probability;
    double *g_PlogM_cum_norm = population.mass_cumulative;

    // ---- Allocate PDMF array [ncomp][nm+1] ----
    double **pdmf = (double**)calloc(ncomp, sizeof(double*));
    for (int i = 0; i < ncomp; i++) pdmf[i] = (double*)calloc(nm+1, sizeof(double));

    // ---- Thin disk (components 0-6): exponential SFR, age bins ----
    {
        int iages7[7] = {15, 100, 200, 300, 500, 700, 1000};
        double gamma  = 1.0 / ctx.stellar.tSFR;
        for (int j = 1; j <= 1000; j++) {
            int itmp = (int)((j - ctx.stellar.agesD[0]) / (double)(ctx.stellar.agesD[1] - ctx.stellar.agesD[0]) + 0.5);
            if (itmp < 0)              itmp = 0;
            if (itmp >= ctx.stellar.nageD) itmp = ctx.stellar.nageD - 1;
            double logMdie = log10(ctx.stellar.MinidieD[itmp]);
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
            double norm = (sum > 0) ? ctx.density.n0MSd[i] / sum : 0;
            for (int iM = 0; iM <= nm; iM++) pdmf[i][iM] *= norm;
        }
    }

    // ---- Thick disk (component 7): single old population ----
    {
        double logMdie = log10(ctx.stellar.MinidieD[ctx.stellar.nageD-2]);
        int iMdie = (int)((logMdie - logMst) / dlogM);
        if (iMdie > nm) iMdie = nm;
        for (int iM = 0; iM <= iMdie && iM <= nm; iM++) pdmf[7][iM] = g_PlogM[iM];
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[7][iM] * dlogM;
        double norm = (sum > 0) ? ctx.density.n0MSd[7] / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[7][iM] *= norm;
    }

    // ---- Bar (component 8): Gaussian SFR around mageB ----
    {
        for (int i = 0; i < ctx.stellar.nageB; i++) {
            double tau  = 0.01 * ctx.stellar.agesB[i];
            double dw   = (tau - ctx.stellar.mageB) / ctx.stellar.sageB;
            double wt   = exp(-0.5 * dw * dw);
            double logMdie = log10(ctx.stellar.MinidieB[i]);
            int iMdie = (int)((logMdie - logMst) / dlogM);
            if (iMdie > nm) iMdie = nm;
            if (iMdie < 0)  continue;
            for (int iM = 0; iM <= iMdie; iM++)
                pdmf[8][iM] += g_PlogM[iM] * wt;
        }
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[8][iM] * dlogM;
        double norm = (sum > 0) ? ctx.density.n0MSb / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[8][iM] *= norm;
    }

    // ---- NSD (component 9): Gaussian SFR around mageND ----
    if (ctx.density.ND > 0 && ctx.density.n0MSND > 0) {
        for (int i = 0; i < ctx.stellar.nageND; i++) {
            double tau  = 0.01 * ctx.stellar.agesND[i];
            double dw   = (tau - ctx.stellar.mageND) / ctx.stellar.sageND;
            double wt   = exp(-0.5 * dw * dw);
            double logMdie = log10(ctx.stellar.MinidieND[i]);
            int iMdie = (int)((logMdie - logMst) / dlogM);
            if (iMdie > nm) iMdie = nm;
            if (iMdie < 0)  continue;
            for (int iM = 0; iM <= iMdie; iM++)
                pdmf[9][iM] += g_PlogM[iM] * wt;
        }
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[9][iM] * dlogM;
        double norm = (sum > 0) ? ctx.density.n0MSND / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[9][iM] *= norm;
    }

    // ---- Stellar halo (component 10): single very old population ----
    if (ctx.density.SH > 0) {
        double logMdie = log10(ctx.stellar.MinidieD[ctx.stellar.nageD-1]);
        int iMdie = (int)((logMdie - logMst) / dlogM);
        if (iMdie > nm) iMdie = nm;
        for (int iM = 0; iM <= iMdie && iM <= nm; iM++) pdmf[10][iM] = g_PlogM[iM];
        double sum = 0;
        for (int iM = 0; iM <= nm; iM++) sum += pdmf[10][iM] * dlogM;
        double norm = (sum > 0) ? ctx.density.n0MSSH / sum : 0;
        for (int iM = 0; iM <= nm; iM++) pdmf[10][iM] *= norm;
    }

    // ---- Header ----
    printf("# calc_mass_dist\n");
    printf("# DISK=%d  hDISK=%d  model=%d  addX=%d  ND=%d  SH=%d\n",
           ctx.density.DISK, ctx.density.hDISK, ctx.density.model, ctx.density.addX,
           ctx.density.ND, ctx.density.SH);
    printf("# IMF: logM in [%.4f, %.4f] with dlogM=%.5f (%d bins)\n",
           logMst, logMst + nm*dlogM, dlogM, nm);

    if (VERBOSE >= 1) {
        const char *names[11] = {
            "disk0","disk1","disk2","disk3","disk4","disk5","disk6",
            "thick","bar","NSD","halo"
        };
        double n0WD[11], n0RG[11], aveMMS[11];
        for (int i = 0; i < 7; i++)  { n0WD[i]=ctx.density.n0d[i]-ctx.density.n0MSd[i]; n0RG[i]=ctx.density.n0RGd[i]; }
        n0WD[7]=ctx.density.n0d[7]-ctx.density.n0MSd[7]; n0RG[7]=ctx.density.n0RGd[7];
        n0WD[8]=ctx.density.n0b -ctx.density.n0MSb;  n0RG[8]=ctx.density.n0RGb;
        n0WD[9]=ctx.density.n0ND-ctx.density.n0MSND; n0RG[9]=ctx.density.n0RGND;
        n0WD[10]=ctx.density.n0SH-ctx.density.n0MSSH; n0RG[10]=ctx.density.n0RGSH;
        double n0MS_arr[11] = {
            ctx.density.n0MSd[0],ctx.density.n0MSd[1],ctx.density.n0MSd[2],ctx.density.n0MSd[3],
            ctx.density.n0MSd[4],ctx.density.n0MSd[5],ctx.density.n0MSd[6],ctx.density.n0MSd[7],
            ctx.density.n0MSb,ctx.density.n0MSND,ctx.density.n0MSSH
        };
        for (int i = 0; i < ncomp; i++) {
            double sumN = 0, sumMN = 0;
            for (int iM = 0; iM <= nm; iM++) {
                double M = pow(10.0, g_logMass[iM]);
                sumN  += pdmf[i][iM] * dlogM;
                sumMN += M * pdmf[i][iM] * dlogM;
            }
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
    population.release_all(ctx);
    free(ctx.density.lDs);
    free(ctx.density.bDs);
    return 0;
}
