#include "genulens/simulation/mass_function.hpp"

#include <cmath>

#include "genulens/simulation/internal/runtime.hpp"

namespace genulens {

void MassFunction::init_from_population(const PopulationRuntime &pop,
                                         const RunContext &ctx) {
    active_state = const_cast<RunContext*>(&ctx);
    log_mass_  = pop.log_mass;
    mass_prob_ = pop.mass_probability;
    mass_cumu_ = pop.mass_cumulative;
    ptiles_    = pop.mass_percentiles;
    nm_     = nm;
    logMst_ = logMst;
    dlogM_  = dlogM;

    M_emps_   = pop.empirical_masses;
    Mag_emps_ = pop.empirical_magnitudes;
    nMLemp_   = pop.empirical_count;
}

double MassFunction::sample_log_mass(double ran, int nbinMmin, int nbinMmax,
                                      double &imp_weight) const {
    double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
    double scaled_ran;
    if (nbinMmax - nbinMmin > 0) {
        double lo = mass_cumu_[nbinMmin];
        double hi = mass_cumu_[nbinMmax];
        imp_weight = hi - lo;
        scaled_ran = ran * imp_weight + lo;
    } else {
        imp_weight = 1.0;
        scaled_ran = ran;
    }
    int inttmp = (int)(scaled_ran * 20);
    int kst = 1;
    for (int itmp = inttmp; itmp > 0; itmp--) {
        kst = ptiles_[itmp] - 1;
        if (kst > 0) break;
    }
    return getcumu2xist(nm_, log_mass_, mass_cumu_, mass_prob_, scaled_ran, kst, 0);
}

MassFunction::RemnantResult MassFunction::evolve(double initial_mass) const {
    void Mini2Mrem(double *pout, double M, int mean);
    double pout[2] = {};
    Mini2Mrem(pout, initial_mass, 0);
    return {pout[0], (int)pout[1]};
}

double MassFunction::cumulative_at(int bin) const {
    return mass_cumu_[bin];
}

double MassFunction::I_magnitude(double mass) const {
    double getx2y(int n, double *x, double *y, double xin);
    if (mass <= 0.010) return 27.060;
    if (mass >= 1.440) return  2.355;
    return getx2y(nMLemp_, M_emps_, Mag_emps_[1], mass);
}

double MassFunction::K_magnitude(double mass) const {
    double getx2y(int n, double *x, double *y, double xin);
    if (mass <= 0.010) return 30.466;
    if (mass >= 1.440) return  1.535;
    return getx2y(nMLemp_, M_emps_, Mag_emps_[4], mass);
}

} // namespace genulens
