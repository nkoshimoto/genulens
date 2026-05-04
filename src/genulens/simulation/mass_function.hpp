#pragma once

#include "genulens/simulation/run_context.hpp"

namespace genulens {

struct PopulationRuntime;

// Wraps IMF and remnant sampling from PopulationRuntime.
// Borrows pointers; the owning PopulationRuntime must outlive this object.
class MassFunction {
public:
    void init_from_population(const PopulationRuntime &pop, const RunContext &ctx);

    // Sample log10(M_initial) from the IMF CDF in bin range [nbinMmin, nbinMmax).
    // When nbinMmin == nbinMmax == 0 the full range is used.
    // Returns the sampled log10(mass) and sets importance_weight to the fraction
    // of the CDF covered by the range (1 when full range).
    double sample_log_mass(double ran, int nbinMmin, int nbinMmax,
                           double &importance_weight) const;

    // Evolve initial mass to present-day mass. Returns remnant type flag:
    //   0 = MS/BD, 1 = WD, 2 = NS, 3 = BH
    struct RemnantResult { double mass; int type; };
    RemnantResult evolve(double initial_mass) const;

    // Apparent magnitude lookup from empirical ML table.
    double I_magnitude(double mass) const;
    double K_magnitude(double mass) const;

    // Mass grid info (needed by EventSampler for mass-range importance sampling)
    int n_bins()        const { return nm_; }
    double log_M_start() const { return logMst_; }
    double d_log_M()     const { return dlogM_; }
    double cumulative_at(int bin) const;

private:
    RunContext *ctx_   = nullptr;
    double *log_mass_  = nullptr;
    double *mass_prob_ = nullptr;
    double *mass_cumu_ = nullptr;
    int    *ptiles_    = nullptr;
    int     nm_        = 0;
    double  logMst_    = 0.0;
    double  dlogM_     = 0.0;

    double  *M_emps_    = nullptr;
    double **Mag_emps_  = nullptr;
    int      nMLemp_    = 0;
};

} // namespace genulens
