/* calc_rho_profile.cpp
 * Output the number density profile along a line of sight as a function of distance.
 * Each Galactic component (thin disk x7, thick disk, bar, NSD, halo) is listed separately.
 *
 * Usage (run from the genulens project root):
 *   calc_rho_profile [options]
 *
 * Key options:
 *   l <deg>        Galactic longitude (default 1.0)
 *   b <deg>        Galactic latitude  (default -3.9)
 *   Dmin <pc>      Minimum distance   (default Dstep)
 *   Dmax <pc>      Maximum distance   (default 16000)
 *   Dstep <pc>     Distance step      (default 100)
 *   VERBOSITY <n>  0: minimal header, 1: full header (default 0)
 *   (plus all galaxy model options from genulens, e.g. model, DISK, addX, ...)
 *
 * Output columns (tab-separated):
 *   D[pc]  nMS_disk0..nMS_disk6  nMS_thick  nMS_bar  nMS_NSD  nMS_halo  nMS_total
 *          n_disk0..n_disk6      n_thick     n_bar    n_NSD    n_halo    n_total
 *
 * where nMS = MS star number density [1/pc^3], n = MS+WD number density [1/pc^3]
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
    printf("Output the number density profile along a line of sight.\n\n");
    printf("Key options:\n");
    printf("  l <deg>        Galactic longitude (default 1.0)\n");
    printf("  b <deg>        Galactic latitude  (default -3.9)\n");
    printf("  Dmin <pc>      Minimum distance   (default Dstep)\n");
    printf("  Dmax <pc>      Maximum distance   (default 16000)\n");
    printf("  Dstep <pc>     Distance step      (default 100)\n");
    printf("  VERBOSITY <n>  Verbosity level    (default 0)\n");
    printf("  -h, --help     Show this help and exit\n\n");
    printf("Source-selection options:\n");
    printf("  SOURCE <n>       1: also output source-weighted rhoD_S columns (default 0)\n");
    printf("  Isrange <lo> <hi>  Apparent I_s range (default 14 21 when AIrc is used)\n");
    printf("  Is <mag>         Single I_s value; uses Is +/- Iserr (default Iserr 0.05)\n");
    printf("  AIrc <mag>       I-band extinction to the red clump\n");
    printf("  IsAIrc <mag>     Alias for AIrc\n");
    printf("  VIsrange <lo> <hi>  Apparent (V-I)_s range\n");
    printf("  VIs <color>      Single (V-I)_s value; uses VIs +/- VIserr (default VIserr 0.05)\n");
    printf("  EVIrc <mag>      E(V-I) reddening to the red clump\n");
    printf("  DMrc <mag>       Red-clump distance modulus (default: Nataf+16 relation)\n");
    printf("  hdust <pc>       Dust scale height (default 164)\n");
    printf("  gammaDs <x>      Source-distance fallback weight when no LF cut is active (default 0.5)\n\n");
    printf("Output columns:\n");
    printf("  D[pc]  nMS[0..10]  nMS_tot  n[0..10]  n_tot");
    printf("  [rhoD_S[0..10]  rhoD_S_tot]\n\n");
    printf("Component indices: 0-6 thin disk, 7 thick disk, 8 bar, 9 NSD, 10 halo.\n");
    printf("Resolution is deterministic and controlled by Dstep; AUTOERR is not used.\n");
}

static int get_numeric_option(int argc, char **argv, const char *name, int argno, double *value)
{
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], name) == 0) {
            if (i + argno >= argc) return 0;
            char *endptr = NULL;
            double v = strtod(argv[i + argno], &endptr);
            if (endptr == argv[i + argno] || *endptr != '\0') return 0;
            *value = v;
            return 1;
        }
    }
    return 0;
}

static int has_option(int argc, char **argv, const char *name)
{
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], name) == 0) return 1;
    return 0;
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

    double lSIMU   = getOptiond(argc, argv, "l",         1,  1.0);
    double bSIMU   = getOptiond(argc, argv, "b",         1, -3.9);

    ctx.density.lDs = (double*)malloc(sizeof(double));
    ctx.density.bDs = (double*)malloc(sizeof(double));
    ctx.density.lDs[0] = lSIMU;
    ctx.density.bDs[0] = bSIMU;

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

    NsdMomentRuntime nsd_moments;
    nsd_moments.initialize_if_enabled(ctx);

    // ---- Tool-specific options ----
    double Dmax    = getOptiond(argc, argv, "Dmax",      1, 16000);
    double Dstep   = getOptiond(argc, argv, "Dstep",     1, 100);
    double Dmin    = getOptiond(argc, argv, "Dmin",      1, Dstep);
    int    VERBOSE = getOptiond(argc, argv, "VERBOSITY", 1, 0);
    int    SOURCE  = getOptioni(argc, argv, "SOURCE",    1, 0);
    if (Dstep <= 0.0) { fprintf(stderr, "Dstep must be > 0\n"); return 1; }
    if (Dmin  <= 0.0) { fprintf(stderr, "Dmin must be > 0\n");  return 1; }
    if (Dmax  <  Dmin) { fprintf(stderr, "Dmax must be >= Dmin\n"); return 1; }

    double Isst  = getOptiond(argc, argv, "Isrange",  1, 14.0);
    double Isen  = getOptiond(argc, argv, "Isrange",  2, 21.0);
    double VIsst = getOptiond(argc, argv, "VIsrange", 1,  0.0);
    double VIsen = getOptiond(argc, argv, "VIsrange", 2,  0.0);
    double tmp;
    if (get_numeric_option(argc, argv, "Ismin", 1, &tmp)) Isst = tmp;
    if (get_numeric_option(argc, argv, "Ismax", 1, &tmp)) Isen = tmp;
    if (get_numeric_option(argc, argv, "Is", 1, &tmp)) {
        double Is2, Iserr = getOptiond(argc, argv, "Iserr", 1, 0.05);
        if (get_numeric_option(argc, argv, "Is", 2, &Is2)) { Isst = tmp; Isen = Is2; }
        else { Isst = tmp - Iserr; Isen = tmp + Iserr; }
    }
    if (get_numeric_option(argc, argv, "VIsmin", 1, &tmp)) VIsst = tmp;
    if (get_numeric_option(argc, argv, "VIsmax", 1, &tmp)) VIsen = tmp;
    if (get_numeric_option(argc, argv, "VIs", 1, &tmp)) {
        double VIs2, VIserr = getOptiond(argc, argv, "VIserr", 1, 0.05);
        if (get_numeric_option(argc, argv, "VIs", 2, &VIs2)) { VIsst = tmp; VIsen = VIs2; }
        else { VIsst = tmp - VIserr; VIsen = tmp + VIserr; }
    }

    double AIrc = getOptiond(argc, argv, "AIrc", 1, 0.0);
    if (AIrc == 0.0) AIrc = getOptiond(argc, argv, "IsAIrc", 1, 0.0);
    double EVIrc   = getOptiond(argc, argv, "EVIrc",   1, 0.0);
    double DMrc    = getOptiond(argc, argv, "DMrc",    1, 0.0);
    double hdust   = getOptiond(argc, argv, "hdust",   1, 164.0);
    double gammaDs = getOptiond(argc, argv, "gammaDs", 1, 0.5);

    if (DMrc == 0.0)
        DMrc = 14.3955 - 0.0239 * lSIMU + 0.0122 * fabs(bSIMU) + 0.128;
    double sinb   = sin(bSIMU / 180.0 * 3.1415926535897932385);
    double hscale = hdust / (fabs(sinb) + 0.0001);
    double Dmean  = (DMrc > 0) ? pow(10.0, 0.2 * DMrc) * 10.0 : -9.99;
    double AI0    = (hscale > 0 && Dmean > 0) ? AIrc  / (1.0 - exp(-Dmean / hscale)) : 0.0;
    double EVI0   = (hscale > 0 && Dmean > 0) ? EVIrc / (1.0 - exp(-Dmean / hscale)) : 0.0;

    int use_I_lf   = (AI0 > 0 && Isen - Isst > 0 && !(EVI0 > 0 && VIsen - VIsst > 0));
    int use_IVI_lf = (AI0 > 0 && EVI0 > 0 && Isen - Isst > 0 && VIsen - VIsst > 0);
    if (use_I_lf || use_IVI_lf || has_option(argc, argv, "AIrc") || has_option(argc, argv, "IsAIrc"))
        SOURCE = 1;

    // Build luminosity functions (stored in ctx.luminosity)
    if (use_I_lf || use_IVI_lf)
        population.initialize_luminosity_functions(ctx, Isst, Isen, VIsst, VIsen, AIrc, EVIrc);

    // ---- Header ----
    printf("# calc_rho_profile\n");
    printf("# (l, b) = (%.4f, %.4f) deg\n", lSIMU, bSIMU);
    printf("# Dmin = %.0f pc, Dmax = %.0f pc, Dstep = %.0f pc\n", Dmin, Dmax, Dstep);
    printf("# DISK=%d  hDISK=%d  model=%d  addX=%d  ND=%d  SH=%d\n",
           ctx.density.DISK, ctx.density.hDISK, ctx.density.model, ctx.density.addX,
           ctx.density.ND, ctx.density.SH);
    if (SOURCE) {
        if (use_IVI_lf)
            printf("# SOURCE=1 with %.3f < Is < %.3f, %.3f < VIs < %.3f\n"
                   "# Extinction: hdust=%.0f pc, DMrc=%.4f, Dmean=%.0f pc, AIrc=%.3f, AI0=%.3f, EVIrc=%.3f, EVI0=%.3f\n",
                   Isst, Isen, VIsst, VIsen, hdust, DMrc, Dmean, AIrc, AI0, EVIrc, EVI0);
        else if (use_I_lf)
            printf("# SOURCE=1 with %.3f < Is < %.3f\n"
                   "# Extinction: hdust=%.0f pc, DMrc=%.4f, Dmean=%.0f pc, AIrc=%.3f, AI0=%.3f\n",
                   Isst, Isen, hdust, DMrc, Dmean, AIrc, AI0);
        else
            printf("# SOURCE=1 without active LF cut; gammaDs=%.3f\n", gammaDs);
    }
    printf("# Columns:\n");
    printf("# D[pc]");
    int ncomp = ctx.density.ncomp;
    for (int i = 0; i < ncomp; i++) printf("  nMS[%d]", i);
    printf("  nMS_tot");
    for (int i = 0; i < ncomp; i++) printf("  n[%d]", i);
    printf("  n_tot");
    if (SOURCE) {
        for (int i = 0; i < ncomp; i++) printf("  rhoD_S[%d]", i);
        printf("  rhoD_S_tot");
    }
    printf("\n");

    auto get_nMS = [&](int i, double shape) -> double {
        if (i < 8)  return ctx.density.n0MSd[i] * shape;
        if (i == 8) return ctx.density.n0MSb    * shape;
        if (i == 9) return ctx.density.n0MSND   * shape;
        return              ctx.density.n0MSSH  * shape;
    };
    auto get_n = [&](int i, double shape) -> double {
        if (i < 8)  return ctx.density.n0d[i] * shape;
        if (i == 8) return ctx.density.n0b    * shape;
        if (i == 9) return ctx.density.n0ND   * shape;
        return              ctx.density.n0SH  * shape;
    };

    double *rhos = (double*)calloc(ncomp, sizeof(double));
    double  xyz[3] = {}, xyb[2] = {};

    for (double D = Dmin; D <= Dmax + 0.5*Dstep; D += Dstep) {
        calc_rho_each(ctx, D, 0, rhos, xyz, xyb);

        double nMS_tot = 0, n_tot = 0;
        printf("%.1f", D);
        for (int i = 0; i < ncomp; i++) {
            double nMS_i = get_nMS(i, rhos[i]);
            nMS_tot += nMS_i;
            printf("\t%.6e", nMS_i);
        }
        printf("\t%.6e", nMS_tot);
        for (int i = 0; i < ncomp; i++) {
            double n_i = get_n(i, rhos[i]);
            n_tot += n_i;
            printf("\t%.6e", n_i);
        }
        printf("\t%.6e", n_tot);
        if (SOURCE) {
            double rhoD_S_tot = 0;
            double extI  = AI0  * (1.0 - exp(-D / hscale)) + 5.0 * log10(0.1 * (D + 0.1));
            double extVI = EVI0 * (1.0 - exp(-D / hscale));
            for (int i = 0; i < ncomp; i++) {
                double nMS_i = get_nMS(i, rhos[i]);
                double rhoD_S = 0;
                if (use_IVI_lf)
                    rhoD_S = nMS_i * fIVI_detect(ctx, extI, Isst, Isen, extVI, VIsst, VIsen, i) * 1e-06 * D * D;
                else if (use_I_lf)
                    rhoD_S = nMS_i * fLF_detect(ctx, extI, Isst, Isen, i) * 1e-06 * D * D;
                else {
                    double w = (gammaDs == 0.5) ? sqrt(D / 8000.0) : pow(((D + 10.0) / 8000.0), fabs(gammaDs));
                    if (gammaDs < 0) w = 1.0 / w;
                    rhoD_S = nMS_i * w * 1e-03;
                }
                rhoD_S_tot += rhoD_S;
                printf("\t%.6e", rhoD_S);
            }
            printf("\t%.6e", rhoD_S_tot);
        }
        printf("\n");
    }

    free(rhos);
    population.release_luminosity_functions(ctx);
    population.release_all(ctx);
    nsd_moments.release_if_enabled(ctx);
    free(ctx.density.lDs);
    free(ctx.density.bDs);
    return 0;
}
