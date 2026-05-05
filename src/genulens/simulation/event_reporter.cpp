#include "genulens/simulation/event_reporter.hpp"

#include <cmath>
#include <cstdio>

namespace genulens {

namespace {

constexpr double STR2MIN2 = 8.461595e-08;

} // namespace

CliEventReporter::CliEventReporter(bool enabled, int verbosity, int binary)
    : enabled_(enabled), verbosity_(verbosity), binary_(binary)
{
}

void CliEventReporter::print_monte_carlo_header() const
{
    if (!enabled_) return;
    if (verbosity_ == 2)
        printf("#        wtj           tE       thetaE          piEN          piEE"
               "   D_S         muSl         muSb iS iL fREM");
    if (verbosity_ == 3)
        printf("#        wtj          M_L   D_L   D_S          t_E      theta_E"
               "         pi_E         pi_EN         pi_EE       mu_rel"
               "        mu_Sl        mu_Sb     I_L     K_L iS iL fREM");
    if (verbosity_ == 4)
        printf("#        wtj          M_L   D_L   D_S          t_E      theta_E"
               "         pi_E         pi_EN         pi_EE       mu_rel"
               "        mu_Sl        mu_Sb     I_L     K_L iS iL fREM"
               "   muhel_N      muhel_E        muhel");
    if (verbosity_ == 5)
        printf("#        wtj          M_L   D_L   D_S          t_E      theta_E"
               "         pi_E         pi_EN         pi_EE       mu_rel"
               "        mu_Sl        mu_Sb     I_L     K_L iS iL fREM"
               "   muhel_N      muhel_E        muhel      R_E");
    if (verbosity_ >= 2 && binary_ == 1)
        printf("     q21         M2         aL     aLpmin         u0 BL");
    if (verbosity_ >= 2)
        printf("\n");
}

void CliEventReporter::print_event(const EventOutputRow &row) const
{
    if (!enabled_ || verbosity_ <= 0) return;
    if (verbosity_ == 1)
        printf("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f"
               " %7.3f %1d %1d %5.2f %5.2f %d %.5e",
               row.wtj, row.tE, row.thetaE, row.piE, row.lens_mass,
               row.source_distance, row.lens_distance,
               row.vx_s, row.vy_s, row.vz_s, row.vx_l, row.vy_l, row.vz_l,
               row.mu_rel, row.source_component, row.lens_component,
               row.tau_s, row.tau_l, row.remnant_flag, row.gamma);
    if (verbosity_ == 2)
        printf("%12.6e %12.6e %12.6e %13.6e %13.6e %5.0f %12.5e %12.5e %2d %2d %d",
               row.wtj, row.tE, row.thetaE, row.piEN, row.piEE,
               row.source_distance, row.muSl, row.muSb,
               row.source_component, row.lens_component, row.remnant_flag);
    if (verbosity_ == 3)
        printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e"
               " %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d",
               row.wtj, row.lens_mass, row.lens_distance, row.source_distance,
               row.tE, row.thetaE, row.piE, row.piEN, row.piEE,
               row.mu_rel, row.muSl, row.muSb, row.IL, row.KL,
               row.source_component, row.lens_component, row.remnant_flag);
    if (verbosity_ == 4)
        printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e"
               " %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d"
               " %12.5e %12.5e %12.6e",
               row.wtj, row.lens_mass, row.lens_distance, row.source_distance,
               row.tE, row.thetaE, row.piE, row.piEN, row.piEE,
               row.mu_rel, row.muSl, row.muSb, row.IL, row.KL,
               row.source_component, row.lens_component, row.remnant_flag,
               row.muhelN, row.muhelE,
               sqrt(row.muhelN*row.muhelN + row.muhelE*row.muhelE));
    if (verbosity_ == 5)
        printf("%12.6e %12.6e %5.0f %5.0f %12.6e %12.6e %12.6e %13.6e %13.6e"
               " %12.6e %12.5e %12.5e %7.3f %7.3f %2d %2d %d"
               " %12.5e %12.5e %12.6e %.6e",
               row.wtj, row.lens_mass, row.lens_distance, row.source_distance,
               row.tE, row.thetaE, row.piE, row.piEN, row.piEE,
               row.mu_rel, row.muSl, row.muSb, row.IL, row.KL,
               row.source_component, row.lens_component, row.remnant_flag,
               row.muhelN, row.muhelE,
               sqrt(row.muhelN*row.muhelN + row.muhelE*row.muhelE),
               row.thetaE_lens_radius);
    if (binary_ == 1)
        printf(" %.4e %.4e %.4e %.4e %.4e %2d",
               row.q21, row.secondary_mass, row.projected_separation,
               row.projected_separation_min, row.u0_source, row.binary_state);
    printf("\n");
}

void CliEventReporter::print_summary(const MonteCarloStats &stats,
                                     long n_simu,
                                     int calc_prior_piE,
                                     int calc_prior_thetaE,
                                     double total_source_count,
                                     double all_source_density,
                                     double tauall) const
{
    if (!enabled_) return;

    printf("# Source number density= %.5e ( %.5e ) arcmin^-2\n",
           total_source_count, all_source_density * STR2MIN2 * 1e+6);

    double medtE = 0;
    double CumuNlogtE = 0;
    for (int ilogtE = 0; ilogtE < stats.NbintE; ilogtE++) {
        double dP = (ilogtE == 0) ? 0.5 * stats.NlogtEs[ilogtE] / stats.SumGamma
                  : 0.5 * (stats.NlogtEs[ilogtE-1] + stats.NlogtEs[ilogtE]) / stats.SumGamma;
        CumuNlogtE += dP;
        if (CumuNlogtE > 0.5) {
            double p2 = CumuNlogtE, p1 = CumuNlogtE - dP;
            double dlogtE = (stats.logtEmax - stats.logtEmin) / stats.NbintE;
            double medlogtE = stats.logtEmin + (ilogtE - 0.5)*dlogtE
                            + (0.5 - p1)/(p2 - p1)*dlogtE;
            medtE = pow(10.0, medlogtE);
            break;
        }
    }
    double avetE       = stats.SumtE / stats.SumGamma;
    double everate     = 2 * tauall / M_PI / avetE * 365.25;
    double everatedeg2 = everate * total_source_count * 3600;
    printf("# avetE= %6.3f days, medtE= %6.3f days, tau= %.6e ,"
           " event_rate= %.6e /star/yr or %.6e /deg^2/yr\n",
           avetE, medtE, tauall, everate, everatedeg2);

    if (binary_ == 1)
        printf("# (n_single n_binwide n_binclose)/n_all="
               " ( %6.0f %6.0f %6.0f ) / %6.0f = ( %.6f %.6f %.6f )\n",
               stats.ncnts, stats.ncntbWD, stats.ncntbCD, stats.ncntall,
               stats.ncnts/stats.ncntall, stats.ncntbWD/stats.ncntall,
               stats.ncntbCD/stats.ncntall);
    printf("# (n_BD n_MS n_WD n_NS n_BH)/n_all="
           " ( %6.0f %6.0f %6.0f %6.0f %6.0f ) / %6.0f"
           " = ( %.6f %.6f %.6f %.6f %.6f )\n",
           stats.nBD, stats.nMS, stats.nWD, stats.nNS, stats.nBH, stats.ncntall,
           stats.nBD/stats.ncntall, stats.nMS/stats.ncntall,
           stats.nWD/stats.ncntall, stats.nNS/stats.ncntall, stats.nBH/stats.ncntall);
    if (stats.NrejIS > 0)
        printf("# N_IS/Ngen= %ld / %ld = %f of events are rejected by importance sampling\n",
               stats.NrejIS, stats.Ngen, (double)stats.NrejIS / stats.Ngen);
    printf("# Nlike/N/Ngen= %d / %ld / %ld     wtlike/wtlike_tE= %.0f / %.0f = %f\n",
           stats.Nlike, n_simu, stats.Ngen,
           stats.wtlike, stats.wtlike_tE,
           stats.wtlike / stats.wtlike_tE);
    if (calc_prior_thetaE)
        printf("# P_pri(thE)= wtlike_w_thEe/wtlike_except_thE= %.2f / %.2f = %.5e\n",
               stats.wtlike_w_thEe, stats.wtlike_except_thE,
               stats.wtlike_w_thEe / stats.wtlike_except_thE);
    if (calc_prior_piE)
        printf("# P_pri(piEN,piEE)= wtlike_w_piEe/wtlike_except_piE= %.2f / %.2f = %.5e\n",
               stats.wtlike_w_piEe, stats.wtlike_except_piE,
               stats.wtlike_w_piEe / stats.wtlike_except_piE);
}

} // namespace genulens
