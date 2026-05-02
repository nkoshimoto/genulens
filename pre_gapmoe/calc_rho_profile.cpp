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

static int nMIs = 0, nVIs = 0;
static double *MIs = NULL, **CumuN_MIs = NULL, dILF = 0;
static double *VIs = NULL, ***f_VI_Is = NULL, dVILF = 0;

static const char *iso_file_for_component(int icomp)
{
    return (icomp == 0) ? "input_files/iso-thin1.dat" :
           (icomp == 1) ? "input_files/iso-thin2.dat" :
           (icomp == 2) ? "input_files/iso-thin3.dat" :
           (icomp == 3) ? "input_files/iso-thin4.dat" :
           (icomp == 4) ? "input_files/iso-thin5.dat" :
           (icomp == 5) ? "input_files/iso-thin6.dat" :
           (icomp == 6) ? "input_files/iso-thin7.dat" :
           (icomp == 7) ? "input_files/iso-thick2.dat" :
           (icomp == 8) ? "input_files/iso-bar_age.dat" :
           (icomp == 9) ? "input_files/iso-NSD.dat" :
                          "input_files/iso-halo.dat";
}

static int iso_nage_for_component(int icomp)
{
    return (icomp == 0) ? 3  : (icomp == 1) ? 18 : (icomp == 2) ? 21 :
           (icomp == 3) ? 21 : (icomp == 4) ? 41 : (icomp == 5) ? 41 :
           (icomp == 6) ? 61 : (icomp == 7) ? 2  : (icomp == 8) ? 27 :
           (icomp == 9) ? 6  : 2;
}

static int iso_narry_for_component(int icomp)
{
    return (icomp == 0) ? 465 : (icomp == 1) ? 581 : (icomp == 2) ? 1885 :
           (icomp == 3) ? 552 : (icomp == 4) ? 362 : (icomp == 5) ? 323 :
           (icomp == 6) ? 296 : (icomp == 7) ? 197 : (icomp == 8) ? 766 :
           (icomp == 9) ? 340 : 190;
}

static int iso_dtau_for_component(int icomp)
{
    return (icomp == 0) ? 10  : (icomp == 1) ? 85  : (icomp == 2) ? 100 :
           (icomp == 3) ? 100 : (icomp == 4) ? 100 : (icomp == 5) ? 50  :
           (icomp == 6) ? 100 : 100;
}

static double component_sfr_weight(int icomp, double tau)
{
    double gamma = 1.0 / tSFR;
    return (icomp == 9) ? exp(-0.5 * pow((tau - mageND) / sageND, 2)) :
           (icomp == 8) ? exp(-0.5 * pow((tau - mageB) / sageB, 2)) :
           (icomp < 7)  ? exp(-gamma * (10 - tau)) :
                          2.0;
}

static int make_LFs(double *MIs_out, double **CumuN_MIs_out, double *PlogM_cum_norm)
{
    FILE *fp;
    char line[1000];
    char *words[100];
    int nMI = 0;
    const int Nbin = 950;
    const int MIst = -6, MIen = 13;
    const double dI = (double)(MIen - MIst) / Nbin;

    if ((fp = fopen("input_files/NbleNall_bin.dat", "r")) == NULL) {
        printf("can't open input_files/NbleNall_bin.dat\n");
        exit(1);
    }
    while (fgets(line, 1000, fp) != NULL) {
        int nwords = split((char*)" ", line, words);
        if (nwords == 0 || *words[0] == '#') continue;
        MIs_out[nMI++] = atof(words[0]);
    }
    fclose(fp);

    for (int icomp = 0; icomp < ncomp; icomp++) {
        const char *file1 = iso_file_for_component(icomp);
        if ((fp = fopen(file1, "r")) == NULL) {
            printf("can't open %s\n", file1);
            exit(1);
        }

        int nage = iso_nage_for_component(icomp);
        int narry = iso_narry_for_component(icomp);
        int dtau = iso_dtau_for_component(icomp);
        int *nMinis = (int*)calloc(nage, sizeof(int));
        double **Minis = (double**)malloc(sizeof(double*) * nage);
        double **MIcs = (double**)malloc(sizeof(double*) * nage);
        for (int j = 0; j < nage; j++) {
            Minis[j] = (double*)calloc(narry, sizeof(double));
            MIcs[j] = (double*)calloc(narry, sizeof(double));
        }

        int iage = 0, iagest = 0;
        while (fgets(line, 1000, fp) != NULL) {
            int nwords = split((char*)" ", line, words);
            if (nwords == 0 || *words[0] == '#') continue;
            if (iagest == 0) iagest = atoi(words[0]);
            if ((atoi(words[0]) - iagest) % dtau != 0) continue;
            iage = (int)((atoi(words[0]) - iagest + 0.5) / dtau);
            if (iage < 0 || iage >= nage) continue;
            Minis[iage][nMinis[iage]] = atof(words[1]);
            MIcs[iage][nMinis[iage]] = atof(words[3]);
            nMinis[iage]++;
        }
        fclose(fp);

        double pIs[960] = {};
        for (int j = 0; j < iage + 1; j++) {
            if (nMinis[j] <= 0) continue;
            double tau = (j * dtau + iagest) * 0.01;
            double wtSFR = component_sfr_weight(icomp, tau);
            if (j == 0 || j == iage) wtSFR *= 0.5;
            if (j == 0 && icomp == 0) wtSFR *= 3.0;
            double PBD = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][0]));
            pIs[Nbin - 5] += wtSFR * PBD;
            for (int k = 0; k < nMinis[j] - 1; k++) {
                if (Minis[j][k + 1] == 0) continue;
                double MIc = 0.5 * (MIcs[j][k] + MIcs[j][k + 1]);
                double P1 = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][k]));
                double P2 = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][k + 1]));
                int pI = (int)((MIc - MIst) / dI);
                if (pI < 0) pI = 0;
                if (pI >= Nbin) pI = Nbin - 5;
                pIs[pI] += wtSFR * (P2 - P1);
            }
        }

        for (int pI = 0; pI <= Nbin; pI++)
            CumuN_MIs_out[icomp][pI] = (pI >= 1) ? 0.5 * (pIs[pI] + pIs[pI - 1]) + CumuN_MIs_out[icomp][pI - 1] : 0.0;

        for (int j = 0; j < nage; j++) {
            free(Minis[j]);
            free(MIcs[j]);
        }
        free(Minis);
        free(MIcs);
        free(nMinis);
    }

    for (int k = 0; k < ncomp; k++)
        for (int j = 0; j < nMI; j++)
            CumuN_MIs_out[k][j] /= CumuN_MIs_out[k][nMI - 1];

    return nMI;
}

static void store_VI_MI(double MIst, double MIen, int NbinMI, double VIst, double VIen, int NbinVI,
                        double *MIs_out, double *VIs_out, double ***f_VI_Is_out, double *PlogM_cum_norm)
{
    FILE *fp;
    char line[1000];
    char *words[100];
    double dMI = (MIen - MIst) / NbinMI;
    double dVI = (VIen - VIst) / NbinVI;

    for (int icomp = 0; icomp < ncomp; icomp++) {
        const char *file1 = iso_file_for_component(icomp);
        if ((fp = fopen(file1, "r")) == NULL) {
            printf("can't open %s\n", file1);
            exit(1);
        }

        int nage = iso_nage_for_component(icomp);
        int narry = iso_narry_for_component(icomp);
        int *nMinis = (int*)calloc(nage, sizeof(int));
        double *ages = (double*)calloc(nage, sizeof(double));
        double **Minis = (double**)malloc(sizeof(double*) * nage);
        double **VIcs = (double**)malloc(sizeof(double*) * nage);
        double **MIcs = (double**)malloc(sizeof(double*) * nage);
        for (int j = 0; j < nage; j++) {
            Minis[j] = (double*)calloc(narry, sizeof(double));
            VIcs[j] = (double*)calloc(narry, sizeof(double));
            MIcs[j] = (double*)calloc(narry, sizeof(double));
        }

        int iage = -1;
        double age0 = -1;
        while (fgets(line, 1000, fp) != NULL) {
            int nwords = split((char*)" ", line, words);
            if (nwords == 0 || *words[0] == '#') continue;
            double age = atof(words[0]);
            if (age != age0) iage++;
            age0 = age;
            if (iage < 0 || iage >= nage) continue;
            Minis[iage][nMinis[iage]] = atof(words[1]);
            VIcs[iage][nMinis[iage]] = atof(words[2]) - atof(words[3]);
            MIcs[iage][nMinis[iage]] = atof(words[3]);
            if (nMinis[iage] == 0) ages[iage] = age;
            nMinis[iage]++;
        }
        fclose(fp);

        double sumwt = 0;
        for (int j = 0; j < iage + 1; j++) {
            if (nMinis[j] <= 0) continue;
            double tau = 0.01 * ages[j];
            double wtSFR = component_sfr_weight(icomp, tau);
            if (j == 0 || j == iage) wtSFR *= 0.5;
            if (j == 0 && icomp == 0) wtSFR *= 3.0;
            double PBD = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][0]));
            f_VI_Is_out[icomp][NbinVI - 5][NbinMI - 5] += wtSFR * PBD;
            sumwt += wtSFR * PBD;
            for (int k = 0; k < nMinis[j] - 1; k++) {
                if (Minis[j][k + 1] == 0) continue;
                double P1 = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][k]));
                double P2 = interp_x(nm + 1, PlogM_cum_norm, logMst, dlogM, log10(Minis[j][k + 1]));
                double wtM = P2 - P1;
                for (int l = 0; l < 3; l++) {
                    double VI = (l == 0) ? VIcs[j][k] : (l == 1) ? 0.5 * (VIcs[j][k] + VIcs[j][k + 1]) : VIcs[j][k + 1];
                    double MI = (l == 0) ? MIcs[j][k] : (l == 1) ? 0.5 * (MIcs[j][k] + MIcs[j][k + 1]) : MIcs[j][k + 1];
                    int pMI = (int)((MI - MIst) / dMI);
                    int pVI = (int)((VI - VIst) / dVI);
                    if (pMI < 0) pMI = 0;
                    if (pMI >= NbinMI) pMI = NbinMI - 1;
                    if (pVI < 0) pVI = 0;
                    if (pVI >= NbinVI) pVI = NbinVI - 1;
                    double tmpwt = (l == 1) ? 1.0 : 0.5;
                    f_VI_Is_out[icomp][pVI][pMI] += tmpwt * wtSFR * wtM;
                    sumwt += tmpwt * wtSFR * wtM;
                }
            }
        }

        for (int pVI = 0; pVI <= NbinVI; pVI++) {
            if (icomp == 0) VIs_out[pVI] = VIst + pVI * dVI;
            for (int pMI = 0; pMI <= NbinMI; pMI++) {
                if (icomp == 0) MIs_out[pMI] = MIst + pMI * dMI;
                f_VI_Is_out[icomp][pVI][pMI] /= sumwt;
            }
        }

        for (int j = 0; j < nage; j++) {
            free(Minis[j]);
            free(VIcs[j]);
            free(MIcs[j]);
        }
        free(Minis);
        free(VIcs);
        free(MIcs);
        free(ages);
        free(nMinis);
    }
}

static double fLF_detect(double extI, double Imin, double Imax, int idisk)
{
    double imaxd = (Imax - extI - MIs[0]) / dILF;
    double imind = (Imin - extI - MIs[0]) / dILF;
    if (imaxd < 0) imaxd = 0;
    if (imaxd > nMIs - 1) imaxd = nMIs - 1;
    if (imind < 0) imind = 0;
    if (imind > nMIs - 1) imind = nMIs - 1;
    int imax = (int)imaxd;
    int imin = (int)imind;
    double fmax = CumuN_MIs[idisk][imax + 1] * (imaxd - imax) + CumuN_MIs[idisk][imax] * (1 - (imaxd - imax));
    double fmin = CumuN_MIs[idisk][imin + 1] * (imind - imin) + CumuN_MIs[idisk][imin] * (1 - (imind - imin));
    return fmax - fmin;
}

static double fIVI_detect(double extI, double Imin, double Imax, double extVI, double VImin, double VImax, int idisk)
{
    double imaxd = (Imax - extI - MIs[0]) / dILF;
    double imind = (Imin - extI - MIs[0]) / dILF;
    double jmaxd = (VImax - extVI - VIs[0]) / dVILF;
    double jmind = (VImin - extVI - VIs[0]) / dVILF;
    if (imaxd < 0) imaxd = 0;
    if (imaxd > nMIs - 1) imaxd = nMIs - 1;
    if (imind < 0) imind = 0;
    if (imind > nMIs - 1) imind = nMIs - 1;
    if (jmaxd < 0) jmaxd = 0;
    if (jmaxd > nVIs - 1) jmaxd = nVIs - 1;
    if (jmind < 0) jmind = 0;
    if (jmind > nVIs - 1) jmind = nVIs - 1;
    double fIVI = 0;
    for (int j = (int)jmind; j <= (int)jmaxd; j++)
        for (int i = (int)imind; i <= (int)imaxd; i++)
            fIVI += f_VI_Is[idisk][j][i];
    return fIVI;
}

int main(int argc, char **argv)
{
    if (has_help_option(argc, argv)) {
        print_usage(argv[0]);
        return 0;
    }

    // ---- Initialize galaxy model ----
    init_galactic_model(argc, argv, 0); // 0: no Shu DF kinematics needed

    // ---- Tool-specific options ----
    double lSIMU   = getOptiond(argc,argv,"l",         1,  1.0);
    double bSIMU   = getOptiond(argc,argv,"b",         1, -3.9);
    double Dmax    = getOptiond(argc,argv,"Dmax",      1, 16000);
    double Dstep   = getOptiond(argc,argv,"Dstep",     1, 100);
    double Dmin    = getOptiond(argc,argv,"Dmin",      1, Dstep);
    int    VERBOSE = getOptiond(argc,argv,"VERBOSITY", 1, 0);
    int    SOURCE  = getOptioni(argc,argv,"SOURCE",    1, 0);
    if (Dstep <= 0.0) {
        fprintf(stderr, "Dstep must be > 0; got %.6g\n", Dstep);
        return 1;
    }
    if (Dmin <= 0.0) {
        fprintf(stderr, "Dmin must be > 0; got %.6g\n", Dmin);
        return 1;
    }
    if (Dmax < Dmin) {
        fprintf(stderr, "Dmax must be >= Dmin; got Dmax=%.6g Dmin=%.6g\n", Dmax, Dmin);
        return 1;
    }

    double Isst  = getOptiond(argc, argv, "Isrange",  1, 14.0);
    double Isen  = getOptiond(argc, argv, "Isrange",  2, 21.0);
    double VIsst = getOptiond(argc, argv, "VIsrange", 1, 0.0);
    double VIsen = getOptiond(argc, argv, "VIsrange", 2, 0.0);
    double tmp;
    if (get_numeric_option(argc, argv, "Ismin", 1, &tmp)) Isst = tmp;
    if (get_numeric_option(argc, argv, "Ismax", 1, &tmp)) Isen = tmp;
    if (get_numeric_option(argc, argv, "Is", 1, &tmp)) {
        double Is2, Iserr = getOptiond(argc, argv, "Iserr", 1, 0.05);
        if (get_numeric_option(argc, argv, "Is", 2, &Is2)) {
            Isst = tmp;
            Isen = Is2;
        } else {
            Isst = tmp - Iserr;
            Isen = tmp + Iserr;
        }
    }
    if (get_numeric_option(argc, argv, "VIsmin", 1, &tmp)) VIsst = tmp;
    if (get_numeric_option(argc, argv, "VIsmax", 1, &tmp)) VIsen = tmp;
    if (get_numeric_option(argc, argv, "VIs", 1, &tmp)) {
        double VIs2, VIserr = getOptiond(argc, argv, "VIserr", 1, 0.05);
        if (get_numeric_option(argc, argv, "VIs", 2, &VIs2)) {
            VIsst = tmp;
            VIsen = VIs2;
        } else {
            VIsst = tmp - VIserr;
            VIsen = tmp + VIserr;
        }
    }

    double AIrc = getOptiond(argc, argv, "AIrc", 1, 0.0);
    if (AIrc == 0.0) AIrc = getOptiond(argc, argv, "IsAIrc", 1, 0.0);
    double EVIrc   = getOptiond(argc, argv, "EVIrc", 1, 0.0);
    double DMrc    = getOptiond(argc, argv, "DMrc", 1, 0.0);
    double hdust   = getOptiond(argc, argv, "hdust", 1, 164.0);
    double gammaDs = getOptiond(argc, argv, "gammaDs", 1, 0.5);

    if (DMrc == 0.0)
        DMrc = 14.3955 - 0.0239 * lSIMU + 0.0122 * fabs(bSIMU) + 0.128;
    double sinb = sin(bSIMU / 180.0 * PI);
    double hscale = hdust / (fabs(sinb) + 0.0001);
    double Dmean = (DMrc > 0) ? pow(10.0, 0.2 * DMrc) * 10.0 : -9.99;
    double AI0 = (hscale > 0 && Dmean > 0) ? AIrc / (1.0 - exp(-Dmean / hscale)) : 0.0;
    double EVI0 = (hscale > 0 && Dmean > 0) ? EVIrc / (1.0 - exp(-Dmean / hscale)) : 0.0;

    int use_I_lf = (AI0 > 0 && Isen - Isst > 0 && !(EVI0 > 0 && VIsen - VIsst > 0));
    int use_IVI_lf = (AI0 > 0 && EVI0 > 0 && Isen - Isst > 0 && VIsen - VIsst > 0);
    if (use_I_lf || use_IVI_lf || has_option(argc, argv, "AIrc") || has_option(argc, argv, "IsAIrc"))
        SOURCE = 1;

    if (use_I_lf) {
        int narry = 960;
        CumuN_MIs = (double**)malloc(sizeof(double*) * ncomp);
        for (int i = 0; i < ncomp; i++)
            CumuN_MIs[i] = (double*)calloc(narry, sizeof(double));
        MIs = (double*)calloc(narry, sizeof(double));
        nMIs = make_LFs(MIs, CumuN_MIs, g_PlogM_cum_norm);
        dILF = (MIs[nMIs - 1] - MIs[0]) / (nMIs - 1);
    }
    if (use_IVI_lf) {
        double MIst = -5.0, MIen = 10.0, VIst = 0.0, VIen = 3.0;
        nMIs = 150;
        nVIs = 30;
        MIs = (double*)calloc(nMIs + 2, sizeof(double));
        VIs = (double*)calloc(nVIs + 2, sizeof(double));
        f_VI_Is = (double***)malloc(sizeof(double**) * ncomp);
        for (int i = 0; i < ncomp; i++) {
            f_VI_Is[i] = (double**)malloc(sizeof(double*) * (nVIs + 2));
            for (int j = 0; j < nVIs + 2; j++)
                f_VI_Is[i][j] = (double*)calloc(nMIs + 2, sizeof(double));
        }
        store_VI_MI(MIst, MIen, nMIs, VIst, VIen, nVIs, MIs, VIs, f_VI_Is, g_PlogM_cum_norm);
        dILF = (MIen - MIst) / nMIs;
        dVILF = (VIen - VIst) / nVIs;
    }

    // lDs[0] / bDs[0] are set inside init_galactic_model from "l" and "b" options,
    // so they are already correct.

    // ---- Header ----
    printf("# calc_rho_profile\n");
    printf("# (l, b) = (%.4f, %.4f) deg\n", lSIMU, bSIMU);
    printf("# Dmin = %.0f pc, Dmax = %.0f pc, Dstep = %.0f pc\n", Dmin, Dmax, Dstep);
    printf("# DISK=%d  hDISK=%d  model=%d  addX=%d  ND=%d  SH=%d\n",
           DISK, hDISK, model, addX, ND, SH);
    if (SOURCE) {
        if (use_IVI_lf) {
            printf("# SOURCE=1 with %.3f < Is < %.3f, %.3f < VIs < %.3f\n",
                   Isst, Isen, VIsst, VIsen);
            printf("# Extinction: hdust=%.0f pc, DMrc=%.4f, Dmean=%.0f pc, AIrc=%.3f, AI0=%.3f, EVIrc=%.3f, EVI0=%.3f\n",
                   hdust, DMrc, Dmean, AIrc, AI0, EVIrc, EVI0);
        } else if (use_I_lf) {
            printf("# SOURCE=1 with %.3f < Is < %.3f\n", Isst, Isen);
            printf("# Extinction: hdust=%.0f pc, DMrc=%.4f, Dmean=%.0f pc, AIrc=%.3f, AI0=%.3f\n",
                   hdust, DMrc, Dmean, AIrc, AI0);
        } else {
            printf("# SOURCE=1 without active LF cut; gammaDs=%.3f\n", gammaDs);
        }
    }
    printf("# Columns:\n");
    printf("# D[pc]");
    for (int i = 0; i < ncomp; i++) printf("  nMS[%d]", i);
    printf("  nMS_tot");
    for (int i = 0; i < ncomp; i++) printf("  n[%d]", i);
    printf("  n_tot");
    if (SOURCE) {
        for (int i = 0; i < ncomp; i++) printf("  rhoD_S[%d]", i);
        printf("  rhoD_S_tot");
    }
    printf("\n");

    // ---- Normalization per component ----
    // disk 0-7: n0MSd[i], n0d[i]
    // bar   8:  n0MSb,    n0b
    // NSD   9:  n0MSND,   n0ND
    // halo 10:  n0MSSH,   n0SH
    auto get_nMS = [](int i, double shape) -> double {
        if (i < 8)  return n0MSd[i] * shape;
        if (i == 8) return n0MSb    * shape;
        if (i == 9) return n0MSND   * shape;
        return              n0MSSH  * shape;
    };
    auto get_n = [](int i, double shape) -> double {
        if (i < 8)  return n0d[i] * shape;
        if (i == 8) return n0b    * shape;
        if (i == 9) return n0ND   * shape;
        return              n0SH  * shape;
    };

    // ---- Loop over distance ----
    double *rhos = (double*)calloc(ncomp, sizeof(double));
    double  xyz[3] = {}, xyb[2] = {};

    for (double D = Dmin; D <= Dmax + 0.5*Dstep; D += Dstep) {
        calc_rho_each(D, 0, rhos, xyz, xyb);

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
            double extI = AI0 * (1.0 - exp(-D / hscale)) + 5.0 * log10(0.1 * (D + 0.1));
            double extVI = EVI0 * (1.0 - exp(-D / hscale));
            for (int i = 0; i < ncomp; i++) {
                double nMS_i = get_nMS(i, rhos[i]);
                double rhoD_S = 0;
                if (use_IVI_lf) {
                    rhoD_S = nMS_i * fIVI_detect(extI, Isst, Isen, extVI, VIsst, VIsen, i) * 1e-06 * D * D;
                } else if (use_I_lf) {
                    rhoD_S = nMS_i * fLF_detect(extI, Isst, Isen, i) * 1e-06 * D * D;
                } else {
                    double tmpDswt = (gammaDs == 0.5) ? sqrt(D / 8000.0) : pow(((D + 10.0) / 8000.0), fabs(gammaDs));
                    if (gammaDs < 0) tmpDswt = 1.0 / tmpDswt;
                    rhoD_S = nMS_i * tmpDswt * 1e-03;
                }
                rhoD_S_tot += rhoD_S;
                printf("\t%.6e", rhoD_S);
            }
            printf("\t%.6e", rhoD_S_tot);
        }
        printf("\n");
    }

    free(rhos);
    if (CumuN_MIs) {
        for (int i = 0; i < ncomp; i++) free(CumuN_MIs[i]);
        free(CumuN_MIs);
    }
    if (f_VI_Is) {
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < nVIs + 2; j++) free(f_VI_Is[i][j]);
            free(f_VI_Is[i]);
        }
        free(f_VI_Is);
    }
    free(MIs);
    free(VIs);
    gsl_rng_free(r_rng);
    return 0;
}
