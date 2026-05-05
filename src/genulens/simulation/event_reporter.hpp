#pragma once

namespace genulens {

struct MonteCarloStats {
    double ncntall = 0, ncnts = 0, ncntbWD = 0, ncntbCD = 0;
    double nBD = 0, nMS = 0, nWD = 0, nNS = 0, nBH = 0;
    int    Nlike = 0;
    long   NrejIS = 0, Ngen = 0;
    double wtlike = 0, wtlike_tE = 0;
    double wtlike_except_piE = 0, wtlike_w_piEe = 0;
    double wtlike_except_thE = 0, wtlike_w_thEe = 0;
    double SumGamma = 0, SumtE = 0;
    double logtEmin = -1, logtEmax = 2;
    int    NbintE = 300;
    double NlogtEs[500] = {};
};

struct EventOutputRow {
    double wtj = 0;
    double tE = 0;
    double thetaE = 0;
    double piE = 0;
    double piEN = 0;
    double piEE = 0;
    double lens_mass = 0;
    double lens_distance = 0;
    double source_distance = 0;
    double vx_s = 0, vy_s = 0, vz_s = 0;
    double vx_l = 0, vy_l = 0, vz_l = 0;
    double mu_rel = 0;
    double muSl = 0, muSb = 0;
    double IL = 0, KL = 0;
    int source_component = 0;
    int lens_component = 0;
    double tau_s = 0, tau_l = 0;
    int remnant_flag = 0;
    double gamma = 0;
    double muhelN = 0, muhelE = 0;
    double thetaE_lens_radius = 0;
    double q21 = 0, secondary_mass = 0, projected_separation = 0;
    double projected_separation_min = 0, u0_source = 0;
    int binary_state = 0;
};

class CliEventReporter {
public:
    CliEventReporter(bool enabled, int verbosity, int binary);

    void print_monte_carlo_header() const;
    void print_event(const EventOutputRow &row) const;
    void print_summary(const MonteCarloStats &stats,
                       long n_simu,
                       int calc_prior_piE,
                       int calc_prior_thetaE,
                       double total_source_count,
                       double all_source_density,
                       double tauall) const;

private:
    bool enabled_;
    int verbosity_;
    int binary_;
};

} // namespace genulens
