#include "genulens/simulation/los_density_grid.hpp"

#include "genulens/model/source_population_prior.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"

namespace genulens {

namespace {

struct SelectionProbabilityCacheStats {
    std::size_t hits = 0;
    std::size_t misses = 0;
};

struct SelectionProbabilityCacheState {
    std::mutex mutex;
    std::unordered_map<std::string, double> cache;
    SelectionProbabilityCacheStats stats;
};

SelectionProbabilityCacheState &selection_probability_cache_state()
{
    static SelectionProbabilityCacheState state;
    return state;
}

double extinction_for_band(const gmodel::BandExtinction &dust, const std::string &band)
{
    if (band == "Vmag") return dust.v_band;
    if (band == "Imag") return dust.i_band;
    if (band == "Jmag_2mass" || band == "Jmag") return dust.j_band;
    if (band == "Hmag_2mass" || band == "Hmag") return dust.h_band;
    if (band == "Ksmag_2mass" || band == "Kmag" || band == "K_L") return dust.k_band;
    if (band == "F087mag") return dust.f087_band;
    if (band == "F146mag") return dust.f146_band;
    if (band == "F213mag") return dust.f213_band;
    return 0.0;
}

std::vector<gmodel::MagnitudeSelection> source_magnitude_selections(
    double distance_pc,
    const model::ExponentialDustExtinction &extinction,
    const LineOfSightDensityGridConfig &cfg)
{
    std::vector<gmodel::MagnitudeSelection> selections;
    const std::size_t n = cfg.source_selection_bands.size();
    selections.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const bool apparent = cfg.source_selection_apparent_magnitudes.empty() ||
                              cfg.source_selection_apparent_magnitudes[i] != 0;
        double offset = 0.0;
        if (apparent) {
            const auto dust = extinction.at_distance(distance_pc);
            offset = extinction.distance_modulus_term(distance_pc) +
                     extinction_for_band(dust, cfg.source_selection_bands[i]);
        }
        selections.push_back({cfg.source_selection_bands[i],
                              cfg.source_selection_min_magnitudes[i],
                              cfg.source_selection_max_magnitudes[i],
                              offset});
    }
    return selections;
}

std::string selection_probability_cache_key(
    int prior_component,
    double distance_pc,
    double min_initial_mass_msun,
    double max_initial_mass_msun,
    const std::vector<gmodel::MagnitudeSelection> &selections,
    const gmodel::ForwardSourceGenerator &generator)
{
    std::ostringstream out;
    out << std::setprecision(17)
        << generator.cache_key()
        << "|component=" << prior_component
        << "|distance_pc=" << distance_pc
        << "|mass=" << min_initial_mass_msun << ',' << max_initial_mass_msun
        << "|selection=";
    for (const auto &selection : selections) {
        out << selection.band << ','
            << selection.min_magnitude << ','
            << selection.max_magnitude << ','
            << selection.magnitude_offset << ';';
    }
    return out.str();
}

double compute_forward_source_selection_probability(
    const RunContext &ctx,
    int prior_component,
    double distance_pc,
    const std::vector<gmodel::MagnitudeSelection> &selections,
    const LineOfSightDensityGridConfig &cfg)
{
    (void)ctx;
    double probability = 0.0;
    double total_weight = 0.0;
    for (const auto &point : gmodel::SourcePopulationPrior::points_for_component(prior_component)) {
        if (!(point.weight > 0.0)) continue;
        gmodel::ForwardSourceQuery query;
        query.component_index = prior_component;
        query.distance_pc = distance_pc;
        query.min_initial_mass_msun = cfg.source_min_initial_mass_msun;
        query.max_initial_mass_msun = cfg.source_max_initial_mass_msun;
        query.use_default_log_age = false;
        query.log_age = point.log_age;
        query.use_default_metallicity = false;
        query.metallicity_mh = point.metallicity_mh;
        query.magnitude_selections = selections;
        probability += point.weight * cfg.forward_source_generator->selection_probability(query);
        total_weight += point.weight;
    }
    return total_weight > 0.0 ? probability / total_weight : 0.0;
}

double forward_source_selection_probability(const RunContext &ctx,
                                            int component,
                                            double distance_pc,
                                            const LineOfSightDensityGridConfig &cfg)
{
    if (cfg.forward_source_generator == nullptr ||
        cfg.source_selection_bands.empty()) {
        return 1.0;
    }

    const int prior_component = (component == 10) ? 7 : component;
    const auto selections = source_magnitude_selections(distance_pc, *cfg.extinction, cfg);
    if (cfg.forward_source_generator->cache_key().rfind("anonymous|", 0) == 0) {
        return compute_forward_source_selection_probability(
            ctx,
            prior_component,
            distance_pc,
            selections,
            cfg);
    }
    const auto key = selection_probability_cache_key(
        prior_component,
        distance_pc,
        cfg.source_min_initial_mass_msun,
        cfg.source_max_initial_mass_msun,
        selections,
        *cfg.forward_source_generator);

    auto &state = selection_probability_cache_state();
    {
        const std::lock_guard<std::mutex> lock(state.mutex);
        const auto found = state.cache.find(key);
        if (found != state.cache.end()) {
            ++state.stats.hits;
            return found->second;
        }
    }

    const double probability = compute_forward_source_selection_probability(
        ctx,
        prior_component,
        distance_pc,
        selections,
        cfg);

    {
        const std::lock_guard<std::mutex> lock(state.mutex);
        if (state.cache.size() > 200000) state.cache.clear();
        state.cache.emplace(key, probability);
        ++state.stats.misses;
    }
    return probability;
}

SelectionProbabilityCacheStats forward_source_selection_probability_cache_stats()
{
    auto &state = selection_probability_cache_state();
    const std::lock_guard<std::mutex> lock(state.mutex);
    return state.stats;
}

} // namespace

static const double a2toSig25BHs[5] = {0,  0.00279,  0.00985,  0.02054,  0.01901};
static const double a1toSig25BHs[5] = {0, -0.02023, -0.07548, -0.16585, -0.15082};
static const double a0toSig25BHs[5] = {1,  1.01342,  1.05382,  1.05898,  0.69889};

void LineOfSightDensityGrid::build(RunContext &ctx,
                                    const LineOfSightDensityGridConfig &cfg) {

    ncomp_ = ctx.density.ncomp;
    const int Dmax = cfg.Dmax;

    nbin_ = (ctx.density.ND > 0 && fabs(cfg.l) < 0.05 && fabs(cfg.b) < 0.05) ? (int)(0.20*Dmax + 0.5)
          : (ctx.density.ND > 0 && fabs(cfg.l) < 0.10 && fabs(cfg.b) < 0.10) ? (int)(0.10*Dmax + 0.5)
          : (ctx.density.ND > 0)                                               ? (int)(0.04*Dmax + 0.5)
          :                                                           (int)(0.01*Dmax + 0.5);
    dD_    = (double)Dmax / nbin_;
    nallS_ = 0.0;
    uses_forward_source_selection_ =
        cfg.forward_source_generator != nullptr && !cfg.source_selection_bands.empty();
    const auto cache_stats_start = forward_source_selection_probability_cache_stats();
    if (uses_forward_source_selection_ && cfg.progress_messages) {
        std::fprintf(stderr,
                     "genulens: building forward-source selection density grid "
                     "(%d distance bins, %d components, %zu band cut%s)...\n",
                     nbin_ + 1,
                     ncomp_,
                     cfg.source_selection_bands.size(),
                     cfg.source_selection_bands.size() == 1 ? "" : "s");
        std::fflush(stderr);
    }

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

    void calc_rho_each(RunContext &ctx, double D, int idata, double *rhos, double *xyz, double *xyb);
    double fLF_detect(RunContext &ctx, double extI, double Imin, double Imax, int idisk);
    double fIVI_detect(RunContext &ctx, double extI, double Imin, double Imax,
                       double extVI, double VImin, double VImax, int idisk);

    double rhos[11] = {}, xyz[3] = {}, xyb[2] = {};
    const int idata = 0;

    int next_progress_percent = 10;
    for (int ibin = 0; ibin <= nbin_; ibin++) {
        D_[ibin] = (double)ibin / nbin_ * Dmax;
        calc_rho_each(ctx, D_[ibin], idata, rhos, xyz, xyb);
        double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
        double z = xyz[2];

        if (cfg.npri > 0 && ibin % cfg.npri == 0)
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
                double sigW0 = (i < 7) ? ctx.kinematics.sigW10d * pow((ctx.kinematics.medtauds.data()[i]+0.01)/10.01, ctx.kinematics.betaW)
                                       : ctx.kinematics.sigW0td;
                double hsigW = (i < 7) ? ctx.kinematics.hsigWt : ctx.kinematics.hsigWT;
                double sigvbs[3] = {};
                if (cfg.BHhb == 1 && i == 8) {
                    void calc_sigvb(RunContext &ctx, double xb, double yb, double zb, double *sigvbs);
                    calc_sigvb(ctx, xyb[0], xyb[1], z, sigvbs);
                }
                double sigW     = (i == 8) ? sigvbs[2] : sigW0 * exp(-(R - ctx.density.R0) / hsigW);
                double sigW2    = sigW * sigW;
                double sigzkick2 = cfg.vkickBH * cfg.vkickBH * PI / 8;
                double sigzadd  = sqrt(sigW2 + sigzkick2);
                double fvBH     = pow(sigzadd / sigW, cfg.betaBH);
                double RhdBH    = (cfg.fixRhdBH == 1) ? cfg.RhdBH0
                                                       : cfg.RhdBH0 * (1 + sigW2 / sigzkick2);
                double fzdBH    = exp((R - ctx.density.R0) / RhdBH) * fvBH;
                double zd0 = (i < 8) ? ctx.density.zd.data()[i] : 235.344943180979;
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

            double nMS = (i == 8)  ? ctx.density.n0MSb  * rhos[8]
                       : (i == 9)  ? ctx.density.n0MSND * rhos[9]
                       : (i == 10) ? ctx.density.n0MSSH * rhos[10]
                       :             ctx.density.n0MSd.data()[i] * rhos[i];
            double rho = (i == 8)  ? ctx.density.n0b   * rhos[8]
                       : (i == 9)  ? ctx.density.n0ND  * rhos[9]
                       : (i == 10) ? ctx.density.n0SH  * rhos[10]
                       :             ctx.density.n0d.data()[i] * rhos[i];

            const double forward_source_prob =
                forward_source_selection_probability(ctx, i, D_[ibin], cfg);
            if (uses_forward_source_selection_) {
                rhoD_S_[i][ibin] = nMS * forward_source_prob * 1e-06 *
                                    D_[ibin] * D_[ibin];
                nallS_ += rhoD_S_[i][ibin] * dD_;
            } else if (cfg.AI0 > 0 && cfg.Isen - cfg.Isst > 0 &&
                cfg.EVI0 > 0 && cfg.VIsen - cfg.VIsst > 0) {
                double fIVIs = fIVI_detect(ctx, extI, cfg.Isst, cfg.Isen,
                                           extVI, cfg.VIsst, cfg.VIsen, i);
                rhoD_S_[i][ibin] = nMS * fIVIs * 1e-06 * D_[ibin] * D_[ibin];
                nallS_ += rhoD_S_[i][ibin] * dD_;
            } else if (cfg.AI0 > 0 && cfg.Isen - cfg.Isst > 0) {
                double fIs = fLF_detect(ctx, extI, cfg.Isst, cfg.Isen, i);
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

            if (cfg.npri > 0 && ibin % cfg.npri == 0) {
                if (cfg.printrhoS) {
                    printf(" %d: %.1e ", i, rhoD_S_[i][ibin]);
                    printf("( %.2e )", cumu_rho_S_[i][ibin]);
                } else {
                    printf(" %d: %.1e ", i, rho);
                    printf("( %.2e )", cumu_rho_L_[i][ibin]);
                }
            }
        }

        if (cfg.npri > 0 && ibin % cfg.npri == 0) {
            printf(" All: %.1e ", rhosum);
            printf(cfg.printrhoS ? "( %.2e )\n" : "( %.2e )\n",
                   cfg.printrhoS ? cumu_rho_all_S_[ibin] : cumu_rho_all_L_[ibin]);
        }

        if (uses_forward_source_selection_ && cfg.progress_messages && nbin_ > 0) {
            const int percent = static_cast<int>(100.0 * ibin / nbin_ + 0.5);
            if (percent >= next_progress_percent && next_progress_percent < 100) {
                std::fprintf(stderr,
                             "genulens: forward-source selection density grid %d%% complete\n",
                             next_progress_percent);
                std::fflush(stderr);
                next_progress_percent += 10;
            }
        }
    }

    if (uses_forward_source_selection_ && cfg.progress_messages) {
        const auto cache_stats_end = forward_source_selection_probability_cache_stats();
        std::fprintf(stderr,
                     "genulens: forward-source selection density grid complete "
                     "(selected source weight %.6e, cache hits %zu, misses %zu)\n",
                     nallS_,
                     cache_stats_end.hits - cache_stats_start.hits,
                     cache_stats_end.misses - cache_stats_start.misses);
        std::fflush(stderr);
    }

    if (cfg.printBHfac) exit(1);

    if (uses_forward_source_selection_ && !(nallS_ > 0.0)) {
        throw std::runtime_error(
            "forward-source selection leaves zero source density; relax the "
            "magnitude cuts, extinction, photometry, or initial-mass range");
    }

    for (int i = 0; i < ncomp_; i++) {
        double norm_S = cumu_rho_S_[i][nbin_];
        double norm_L = cumu_rho_L_[i][nbin_];
        if (norm_S == 0) continue;
        for (int ibin = 0; ibin <= nbin_; ibin++) {
            double Pnorm_S = cumu_rho_S_[i][ibin] / norm_S;
            int intp_S = (int)(Pnorm_S * 20);
            if (ibinptiles_S_[i][intp_S] == 0)
                ibinptiles_S_[i][intp_S] = (intp_S == 0) ? 1 : (int)(ibin + 0.5);
            if (norm_L == 0) continue;
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
