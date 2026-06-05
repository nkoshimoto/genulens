#include "genulens/simulation/event_sampler.hpp"

#include "genulens/model/source_population_prior.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "genulens/simulation/event_reporter.hpp"
#include "genulens/simulation/internal/runtime.hpp"

#include <limits>
#define MINSIGLOGA 0.3
#define MAXMEANLOGA 1.7
#define MINMEANLOGA 0.6

namespace genulens {
namespace {

bool has_magnitude(const model::ForwardSource &source, const std::string &band)
{
    return source.stellar.absolute_magnitudes.find(band) != source.stellar.absolute_magnitudes.end();
}

double magnitude_or_nan(const model::ForwardSource &source, const std::string &band)
{
    const auto found = source.stellar.absolute_magnitudes.find(band);
    return found == source.stellar.absolute_magnitudes.end()
               ? std::numeric_limits<double>::quiet_NaN()
               : found->second;
}

double extinction_for_band(const model::BandExtinction &dust,
                           const std::string &band)
{
    if (band == "Vmag") return dust.v_band;
    if (band == "Imag") return dust.i_band;
    if (band == "Jmag_2mass" || band == "Jmag") return dust.j_band;
    if (band == "Hmag_2mass" || band == "Hmag") return dust.h_band;
    if (band == "Ksmag_2mass" || band == "Kmag" || band == "K_L") return dust.k_band;
    if (band == "F087mag") return dust.f087_band;
    if (band == "F146mag") return dust.f146_band;
    if (band == "F213mag") return dust.f213_band;
    throw std::runtime_error("apparent magnitude selection has no extinction coefficient for band: " +
                             band);
}

double selected_magnitude(const model::ForwardSource &source,
                          const model::ExponentialDustExtinction *extinction,
                          const std::string &band,
                          bool apparent)
{
    double mag = magnitude_or_nan(source, band);
    if (apparent) {
        if (extinction == nullptr) {
            throw std::runtime_error("apparent magnitude selection requires an extinction model");
        }
        const auto dust = extinction->at_distance(source.distance_pc);
        mag += extinction->distance_modulus_term(source.distance_pc);
        mag += extinction_for_band(dust, band);
    }
    return mag;
}

std::vector<model::MagnitudeSelection> source_magnitude_selections(
    double distance_pc,
    const model::ExponentialDustExtinction *extinction,
    const EventSampler::Config &cfg)
{
    std::vector<model::MagnitudeSelection> selections;
    const std::size_t n = cfg.source_selection_bands.size();
    selections.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const bool apparent = cfg.source_selection_apparent_magnitudes.empty() ||
                              cfg.source_selection_apparent_magnitudes[i] != 0;
        double offset = 0.0;
        if (apparent) {
            if (extinction == nullptr) {
                throw std::runtime_error("apparent magnitude selection requires an extinction model");
            }
            const auto dust = extinction->at_distance(distance_pc);
            offset = extinction->distance_modulus_term(distance_pc) +
                     extinction_for_band(dust, cfg.source_selection_bands[i]);
        }
        selections.push_back({cfg.source_selection_bands[i],
                              cfg.source_selection_min_magnitudes[i],
                              cfg.source_selection_max_magnitudes[i],
                              offset});
    }
    return selections;
}

double source_selection_proposal_distance(double distance_pc,
                                          const EventSampler::Config &cfg)
{
    if (!(cfg.source_selection_distance_bin_pc > 0.0)) return distance_pc;
    const double bin = cfg.source_selection_distance_bin_pc;
    return std::max(0.0, std::round(distance_pc / bin) * bin);
}

std::vector<model::MagnitudeSelection> source_magnitude_proposal_selections(
    double distance_pc,
    const model::ExponentialDustExtinction *extinction,
    const EventSampler::Config &cfg)
{
    if (!(cfg.source_selection_distance_bin_pc > 0.0)) {
        return source_magnitude_selections(distance_pc, extinction, cfg);
    }

    const double bin = cfg.source_selection_distance_bin_pc;
    const double center = source_selection_proposal_distance(distance_pc, cfg);
    const double lo_distance = std::max(1.0e-6, center - 0.5 * bin);
    const double hi_distance = std::max(lo_distance, center + 0.5 * bin);

    std::vector<model::MagnitudeSelection> selections;
    const std::size_t n = cfg.source_selection_bands.size();
    selections.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const bool apparent = cfg.source_selection_apparent_magnitudes.empty() ||
                              cfg.source_selection_apparent_magnitudes[i] != 0;
        double min_absolute = cfg.source_selection_min_magnitudes[i];
        double max_absolute = cfg.source_selection_max_magnitudes[i];
        if (apparent) {
            if (extinction == nullptr) {
                throw std::runtime_error("apparent magnitude selection requires an extinction model");
            }
            const auto lo_dust = extinction->at_distance(lo_distance);
            const auto hi_dust = extinction->at_distance(hi_distance);
            const double lo_offset = extinction->distance_modulus_term(lo_distance) +
                                     extinction_for_band(lo_dust, cfg.source_selection_bands[i]);
            const double hi_offset = extinction->distance_modulus_term(hi_distance) +
                                     extinction_for_band(hi_dust, cfg.source_selection_bands[i]);
            const double min_offset = std::min(lo_offset, hi_offset);
            const double max_offset = std::max(lo_offset, hi_offset);
            min_absolute = cfg.source_selection_min_magnitudes[i] - max_offset;
            max_absolute = cfg.source_selection_max_magnitudes[i] - min_offset;
        }
        selections.push_back({cfg.source_selection_bands[i],
                              min_absolute,
                              max_absolute,
                              0.0});
    }
    return selections;
}

bool passes_source_selection(const model::ForwardSource &source,
                             const model::ExponentialDustExtinction *extinction,
                             const EventSampler::Config &cfg)
{
    if (!cfg.source_selection_bands.empty()) {
        const std::size_t n = cfg.source_selection_bands.size();
        if (cfg.source_selection_min_magnitudes.size() != n ||
            cfg.source_selection_max_magnitudes.size() != n) {
            throw std::runtime_error("forward source selection band/range length mismatch");
        }
        if (!cfg.source_selection_apparent_magnitudes.empty() &&
            cfg.source_selection_apparent_magnitudes.size() != n) {
            throw std::runtime_error("forward source apparent-selection flag length mismatch");
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!has_magnitude(source, cfg.source_selection_bands[i])) return false;
            const bool apparent = cfg.source_selection_apparent_magnitudes.empty() ||
                                  cfg.source_selection_apparent_magnitudes[i] != 0;
            const double mag = selected_magnitude(source, extinction,
                                                  cfg.source_selection_bands[i],
                                                  apparent);
            if (!(mag >= cfg.source_selection_min_magnitudes[i] &&
                  mag <= cfg.source_selection_max_magnitudes[i])) {
                return false;
            }
        }
        return true;
    }

    return true;
}

model::ForwardSource sample_matched_source(const model::ForwardSourceGenerator &generator,
                                           const model::ForwardSourceQuery &query,
                                           const EventSampler::Config &cfg)
{
    auto selected_query = query;
    if (!cfg.source_selection_bands.empty()) {
        selected_query.magnitude_selections =
            source_magnitude_proposal_selections(query.distance_pc, cfg.extinction, cfg);
    }
    const int max_attempts = cfg.source_selection_bands.empty() ? 1 : 16;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        const auto source = generator.sample(selected_query, *cfg.forward_source_rng);
        if (passes_source_selection(source, cfg.extinction, cfg)) {
            return source;
        }
    }
    throw std::runtime_error("forward source did not pass exact source selection");
}

model::AgeMetallicityPoint sample_source_population_point(
    const model::ForwardSourceGenerator &generator,
    int component,
    double distance_pc,
    const EventSampler::Config &cfg)
{
    const int prior_component = (component == 10) ? 7 : component;
    const auto points = model::SourcePopulationPrior::points_for_component(prior_component);
    std::vector<double> cumulative;
    cumulative.reserve(points.size());
    double total = 0.0;
    const auto selections = source_magnitude_proposal_selections(distance_pc, cfg.extinction, cfg);
    for (const auto &point : points) {
        if (!(point.weight > 0.0)) {
            cumulative.push_back(total);
            continue;
        }
        double weight = point.weight;
        if (!selections.empty()) {
            model::ForwardSourceQuery query;
            query.component_index = prior_component;
            query.distance_pc = distance_pc;
            query.min_initial_mass_msun = cfg.source_min_initial_mass_msun;
            query.max_initial_mass_msun = cfg.source_max_initial_mass_msun;
            query.use_default_log_age = false;
            query.log_age = point.log_age;
            query.use_default_metallicity = false;
            query.metallicity_mh = point.metallicity_mh;
            query.magnitude_selections = selections;
            weight *= generator.selection_probability(query);
        }
        total += weight;
        cumulative.push_back(total);
    }
    if (!(total > 0.0)) {
        throw std::runtime_error("forward-source population prior has zero selected probability");
    }
    const double draw = cfg.forward_source_rng->uniform() * total;
    for (std::size_t i = 0; i < cumulative.size(); ++i) {
        if (draw <= cumulative[i]) return points[i];
    }
    return points.back();
}

} // namespace

int EventSampler::run(RunContext &ctx,
	                           LineOfSightDensityGrid &grid,
	                           PopulationRuntime &pop,
	                           MassFunction &mf,
	                           const Config &cfg,
	                           const ObservationConfig &obs,
	                           LikelihoodFunction custom_likelihood,
	                           EventSink event_sink,
	                           bool emit_cli_output) {

    // Sampling options aliases (mutable, mirror RunContext fields)
    auto &samp = ctx.sampling;
    long &NSIMU      = samp.n_simu;
    long &NlikeMIN   = samp.n_like_min;
    int  &VERBOSITY  = samp.verbosity;
    int  &UNIFORM    = samp.uniform_likelihood;
    int  &REMNANT    = samp.remnant;
    int  &BINARY     = samp.binary;
    int  &onlyWD     = samp.only_white_dwarf;
    int  &SMALLGAMMA = samp.small_gamma;
    int  &NoGAMMAIS  = samp.no_gamma_importance_sampling;
    double &wtD_L    = samp.weight_lens_distance;
    double &wtM_L    = samp.weight_lens_mass;
    double &vEarthl  = samp.v_earth_l;
    double &vEarthb  = samp.v_earth_b;
    int  &CALCPRIORpiE = samp.calc_prior_piE;
    int  &CALCPRIORthE = samp.calc_prior_thetaE;

    // Observation shortcuts
    const double tEobs      = obs.tE_obs,    tEe      = obs.tE_err,    fetE     = obs.fe_tE;
    const int    tEdet      = obs.tE_det;
    const double tEmin      = obs.tE_min,    tEmax    = obs.tE_max;
    const double thetaEobs  = obs.thetaE_obs, thetaEe = obs.thetaE_err, fethetaE = obs.fe_thetaE;
    const int    thetaEdet  = obs.thetaE_det;
    const double thetaEmin  = obs.thetaE_min, thetaEmax = obs.thetaE_max;
    const double piEobs     = obs.piE_obs,   piEe     = obs.piE_err,   fepiE    = obs.fe_piE;
    const int    piEdet     = obs.piE_det;
    const double piEmin     = obs.piE_min,   piEmax   = obs.piE_max;
    const double piENobs    = obs.piEN_obs,  piENe    = obs.piEN_err,  fepiEN   = obs.fe_piEN;
    const double piEEobs    = obs.piEE_obs,  piEEe    = obs.piEE_err,  fepiEE   = obs.fe_piEE;
    const double muslobs    = obs.musl_obs,  musle    = obs.musl_err,  femusl   = obs.fe_musl;
    const double musbobs    = obs.musb_obs,  musbe    = obs.musb_err,  femusb   = obs.fe_musb;
    const double musNobs    = obs.musN_obs,  musNe    = obs.musN_err,  femusN   = obs.fe_musN;
    const double musEobs    = obs.musE_obs,  musEe    = obs.musE_err,  femusE   = obs.fe_musE;
    const int    musRCG     = obs.musRCG;
    const double muhelNobs  = obs.muhelN_obs, muhelNe  = obs.muhelN_err, femuhelN = obs.fe_muhelN;
    const double muhelEobs  = obs.muhelE_obs, muhelEe  = obs.muhelE_err, femuhelE = obs.fe_muhelE;
    const double ILobs      = obs.IL_obs,   ILe      = obs.IL_err,   feIL     = obs.fe_IL;
    const int    ILdet      = obs.IL_det;
    const double KLobs      = obs.KL_obs,   KLe      = obs.KL_err,   feKL     = obs.fe_KL;
    const int    KLdet      = obs.KL_det;
    const double u0obs      = obs.u0_obs;

    // Geometry
    const double cosPA = cfg.cosPA, sinPA = cfg.sinPA;
    const double cosb  = cfg.cosb,  sinb  = cfg.sinb;
    const double cosl  = cfg.cosl,  sinl  = cfg.sinl;

    // Extinction
    const auto  *ext  = cfg.extinction;
    const double AI0  = cfg.AI0;
    const double AK0  = cfg.AK0;

    // BH kick
    const int    BHhd    = cfg.BHhd,    BHhb    = cfg.BHhb;
    const double vkickBH = cfg.vkickBH, vkickNS = cfg.vkickNS;
    const int    MXDkick = cfg.MXDkick;
    const double betaBH  = cfg.betaBH;
    const double Dmean   = cfg.Dmean;
    const bool attach_source_properties =
        cfg.forward_source_generator != nullptr && cfg.forward_source_rng != nullptr;

    // PopulationRuntime raw pointers (still used directly)
    double *logMass_B        = pop.log_mass;
    double *PlogM_cum_norm_B = pop.mass_cumulative;
    double *PlogM_B          = pop.mass_probability;
    int    *imptiles_B       = pop.mass_percentiles;

    // Function declarations from runtime
    void get_vxyz_ran(RunContext &ctx, double *vxyz, int i, double tau, double D, double lD, double bD);
    void Mini2Mrem(RunContext &ctx, double *pout, double M, int mean);
    void getaproj(RunContext &ctx, double *pout, double M1, double M2, int coeff);
    double like_obs(RunContext &ctx, double mod, double obs_val, double err, double fe, int det, int uniform_mode);
    double getcumu2xist(int n, double *x, double *F, double *f, double Freq, int ist, int inv);
    double getx2y(int n, double *x, double *y, double xin);
    double interp_x(int n, double *F, double xst, double dx, double xreq);

    CliEventReporter reporter(emit_cli_output, VERBOSITY, BINARY);
    reporter.print_monte_carlo_header();

    MonteCarloStats stats;

    for (long j = 0; j < NSIMU; j++) {
        if (j == NSIMU-1 && NlikeMIN > 0) {
            long NSIMUnew = (NSIMU-1) * (double)NlikeMIN / (stats.Nlike + 1);
	            if (NSIMUnew > NSIMU) { if (emit_cli_output) printf("# Increase NSIMU to %ld\n", NSIMUnew); NSIMU = NSIMUnew; }
	            if (NSIMU > NDATAMAX) { if (emit_cli_output) printf("# Decrease NSIMU to %ld\n", NDATAMAX); NSIMU = NDATAMAX; }
        }

        double ran, cumu, addGammaIS = 1;
        int inttmp, kst;
        stats.Ngen++;

        // --- Sample source component and distance ---
        int i_s = grid.sample_source_component(ctx.runtime.rng->uniform());
        if (i_s == ctx.density.ncomp) { j--; continue; }

        double tau_s = (i_s == 9)  ? ctx.stellar.mageND + ctx.stellar.sageND*ctx.runtime.rng->gaussian()
                     : (i_s == 8)  ? ctx.stellar.mageB  + ctx.stellar.sageB *ctx.runtime.rng->gaussian()
                     : (i_s == 10) ? 14.0
                     :               ctx.kinematics.medtauds.data()[i_s];
        double D_s = grid.sample_source_distance(i_s, ctx.runtime.rng->uniform());
        int nbinDs = (int)(D_s / grid.dD());

        // --- Importance sampling lens distance range ---
        int nbinDlmin = 0, nbinDlmax = nbinDs;
        grid.get_lens_importance_bins(D_s,
                                       thetaEmin, thetaEmax,
                                       piEmin, piEmax,
                                       nbinDlmin, nbinDlmax);

        // --- Sample lens component and distance ---
        int i_l = grid.sample_lens_component(nbinDlmin, nbinDlmax, ctx.runtime.rng->uniform());
        if (i_l == ctx.density.ncomp) { j--; continue; }

        double tau_l = (i_l == 9)  ? ctx.stellar.mageND + ctx.stellar.sageND*ctx.runtime.rng->gaussian()
                     : (i_l == 8)  ? ctx.stellar.mageB  + ctx.stellar.sageB *ctx.runtime.rng->gaussian()
                     : (i_l == 10) ? 14.0
                     :               ctx.kinematics.medtauds.data()[i_l];
        double D_l = grid.sample_lens_distance(i_l, nbinDlmin, nbinDlmax, ctx.runtime.rng->uniform());
        addGammaIS *= grid.importance_weight(nbinDlmin, nbinDlmax);

        // --- Velocities ---
        double vxyz_S[3] = {}, vxyz_L[3] = {};
        get_vxyz_ran(ctx, vxyz_S, i_s, tau_s, D_s, cfg.l, cfg.b);
        get_vxyz_ran(ctx, vxyz_L, i_l, tau_l, D_l, cfg.l, cfg.b);
        double vx_s = vxyz_S[0], vx_l = vxyz_L[0];
        double vy_s = vxyz_S[1], vy_l = vxyz_L[1];
        double vz_s = vxyz_S[2], vz_l = vxyz_L[2];

        double vxrel_s = vx_s - ctx.kinematics.vxsun;
        double vyrel_s = vy_s - ctx.kinematics.vysun;
        double vzrel_s = vz_s - ctx.kinematics.vzsun;
        double muSl = (vxrel_s*sinl      + vyrel_s*cosl)*KS2MY/D_s;
        double muSb = (vxrel_s*cosl*sinb - vyrel_s*sinl*sinb + vzrel_s*cosb)*KS2MY/D_s;

        double murells[3] = {}, murelbs[3] = {}, murels[3] = {};
        double thetakick = (REMNANT == 1 && onlyWD == 0) ? asin(1 - 2*ctx.runtime.rng->uniform()) : 0;
        double phikick   = (REMNANT == 1 && onlyWD == 0) ? ctx.runtime.rng->uniform()*2*PI - PI : 0;
        int nREM = (REMNANT == 1 && onlyWD == 0) ? 3 : 1;
        for (int iREM = 0; iREM < nREM; iREM++) {
            double vxadd = 0, vyadd = 0, vzadd = 0;
            if (iREM > 0) {
                double vkick = (iREM == 1) ? vkickNS : vkickBH;
                if (MXDkick == 1) {
                    double sigv1D = vkick * 0.5 * sqrt(0.5*PI);
                    vxadd = sigv1D * ctx.runtime.rng->gaussian();
                    vyadd = sigv1D * ctx.runtime.rng->gaussian();
                    vzadd = sigv1D * ctx.runtime.rng->gaussian();
                } else {
                    vxadd = vkick * cos(thetakick) * cos(phikick);
                    vyadd = vkick * cos(thetakick) * sin(phikick);
                    vzadd = vkick * sin(thetakick);
                }
            }
            double vxrel_l = vx_l + vxadd - ctx.kinematics.vxsun;
            double vyrel_l = vy_l + vyadd - ctx.kinematics.vysun;
            double vzrel_l = vz_l + vzadd - ctx.kinematics.vzsun;
            double muLl = (vxrel_l*sinl      + vyrel_l*cosl)*KS2MY/D_l;
            double muLb = (vxrel_l*cosl*sinb - vyrel_l*sinl*sinb + vzrel_l*cosb)*KS2MY/D_l;
            double murellhel = muLl - muSl;
            double murelbhel = muLb - muSb;
            murells[iREM] = murellhel - vEarthl*KS2MY*(D_s-D_l)/D_s/D_l;
            murelbs[iREM] = murelbhel - vEarthb*KS2MY*(D_s-D_l)/D_s/D_l;
            murels[iREM]  = sqrt(murells[iREM]*murells[iREM] + murelbs[iREM]*murelbs[iREM]);
        }
        double murell = murells[0];
        double murelb = murelbs[0];
        double murel  = murels[0];

        // --- Mass range from importance sampling ---
        double pirel = 1000.0*(1.0/D_l - 1.0/D_s);

        double MBHmax = 15.6775, MWDmax = 1.375;
        double Minidie;
        if (i_l == 8) {
            int iage_l = (int)(tau_l * 2 + 0.5) * 50;
            int itmp = (iage_l - ctx.stellar.agesB.data()[0]) / (ctx.stellar.agesB.data()[1] - ctx.stellar.agesB.data()[0]);
            Minidie = ctx.stellar.MinidieB.data()[itmp];
        } else if (i_l == 10) {
            Minidie = ctx.stellar.MinidieD.data()[ctx.stellar.nageD-1];
        } else if (i_l == 9) {
            int iage_l = (int)(100*(tau_l + 0.5));
            int itmp = (iage_l <= ctx.stellar.agesND.data()[0]) ? 0
                     : (iage_l >= ctx.stellar.agesND.data()[ctx.stellar.nageND-1]) ? ctx.stellar.nageND-1
                     : (iage_l - ctx.stellar.agesND.data()[0]) / (ctx.stellar.agesND.data()[1] - ctx.stellar.agesND.data()[0]);
            Minidie = ctx.stellar.MinidieND.data()[itmp];
        } else if (i_l == 7) {
            Minidie = ctx.stellar.MinidieD.data()[ctx.stellar.nageD-2];
        } else {
            int iage_l = (int)(tau_l * 100 + 0.5);
            iage_l = (iage_l % 5 > 2) ? iage_l + (5 - iage_l % 5) : iage_l - iage_l % 5;
            if (iage_l < 5) iage_l = 5;
            int itmp = (iage_l - ctx.stellar.agesD.data()[0]) / (ctx.stellar.agesD.data()[1] - ctx.stellar.agesD.data()[0]);
            Minidie = ctx.stellar.MinidieD.data()[itmp];
        }
        double MMSmax = (BINARY == 1) ? 2 * Minidie : Minidie;
        double MPDmax = (REMNANT == 1) ? MBHmax
                       : (onlyWD == 1 && MWDmax > MMSmax) ? MWDmax : MMSmax;
        double MPDmin = ctx.imf_options.ml;

        double MtEmin = MPDmax, MtEmax = MPDmin;
        double MthEmin = MPDmin, MthEmax = MPDmax;
        double MpiEmin = MPDmin, MpiEmax = MPDmax;

        if (tEmax - tEmin > 0) {
            double MMSWDmax = (MWDmax > MMSmax && (onlyWD == 1 || REMNANT == 1)) ? MWDmax : MMSmax;
            double tmpfac = murel*murel / KAPPA/pirel/365.25/365.25;
            double MtEMSmin = tmpfac*tEmin*tEmin; if (MtEMSmin < MPDmin) MtEMSmin = MPDmin;
            double MtEMSmax = tmpfac*tEmax*tEmax; if (MtEMSmax > MMSWDmax) MtEMSmax = MMSWDmax;
            if (MtEMSmax > MPDmin && MtEMSmin < MMSWDmax) { MtEmin = MtEMSmin; MtEmax = MtEMSmax; }
            if (REMNANT == 1 && onlyWD == 0) {
                tmpfac = murels[1]*murels[1] / KAPPA/pirel/365.25/365.25;
                double MtENSmin = tmpfac*tEmin*tEmin; if (MtENSmin < MNSMIN) MtENSmin = MNSMIN;
                double MtENSmax = tmpfac*tEmax*tEmax; if (MtENSmax > MNSMAX) MtENSmax = MNSMAX;
                if (MtENSmax > MNSMIN && MtENSmin < MNSMAX) {
                    if (MtENSmin < MtEmin) MtEmin = MtENSmin;
                    if (MtENSmax > MtEmax) MtEmax = MtENSmax;
                }
                double MBHmin = 4.9911;
                tmpfac = murels[2]*murels[2] / KAPPA/pirel/365.25/365.25;
                double MtEBHmin = tmpfac*tEmin*tEmin; if (MtEBHmin < MBHmin) MtEBHmin = MBHmin;
                double MtEBHmax = tmpfac*tEmax*tEmax; if (MtEBHmax > MBHmax) MtEBHmax = MBHmax;
                if (MtEBHmax > MBHmin && MtEBHmin < MBHmax) {
                    if (MtEBHmin < MtEmin) MtEmin = MtEBHmin;
                    if (MtEBHmax > MtEmax) MtEmax = MtEBHmax;
                }
            }
            if (MtEmin > MtEmax) { j--; stats.NrejIS++; continue; }
        } else {
            MtEmin = MPDmin; MtEmax = MPDmax;
        }

        if (thetaEmax - thetaEmin > 0) {
            double Mtmp;
            Mtmp = thetaEmin*thetaEmin/KAPPA/pirel;
            MthEmin = (Mtmp < MPDmin) ? MPDmin : (Mtmp > MPDmax) ? MPDmax : Mtmp;
            Mtmp = thetaEmax*thetaEmax/KAPPA/pirel;
            MthEmax = (Mtmp < MPDmin) ? MPDmin : (Mtmp > MPDmax) ? MPDmax : Mtmp;
        }
        if (piEmax - piEmin > 0) {
            double Mtmp;
            Mtmp = pirel/KAPPA/piEmax/piEmax;
            MpiEmin = (Mtmp < MPDmin) ? MPDmin : (Mtmp > MPDmax) ? MPDmax : Mtmp;
            Mtmp = pirel/KAPPA/piEmin/piEmin;
            MpiEmax = (Mtmp < MPDmin) ? MPDmin : (Mtmp > MPDmax) ? MPDmax : Mtmp;
        }

        if (MtEmin >= MthEmax || MtEmin >= MpiEmax || MthEmin >= MpiEmax ||
            MtEmax <= MthEmin || MtEmax <= MpiEmin || MthEmax <= MpiEmin ||
            MtEmax - MtEmin <= 0 || MthEmax - MthEmin <= 0 || MpiEmax - MpiEmin <= 0) {
            j--; stats.NrejIS++; continue;
        }

        // --- Sample lens mass ---
        int nbinMmin = 0, nbinMmax = 0;
        if (tEmax - tEmin > 0 || thetaEmax - thetaEmin > 0 || piEmax - piEmin > 0) {
            double Minimin = (MtEmin > MthEmin) ? MtEmin : MthEmin;
            if (MpiEmin > Minimin) Minimin = MpiEmin;
            double Minimax = (MtEmax < MthEmax) ? MtEmax : MthEmax;
            if (MpiEmax < Minimax) Minimax = MpiEmax;
            if (REMNANT == 1 && onlyWD == 0) {
                double MWDmin = 0.109 * Minidie + 0.394;
                Minimax = (Minimax > MNSMIN) ? ctx.imf_options.mu
                        : (Minimax > MWDmin) ? (Minimax - 0.394)/0.109
                        : Minimax;
                Minimin = (Minimin > MWDmax && Minimin > MMSmax) ? ctx.stellar.MiniWDmax
                        : (Minimin > MMSmax) ? (Minimin - 0.394)/0.109
                        : (BINARY == 1) ? 0.5 * Minimin : Minimin;
            } else if (onlyWD == 1) {
                double MWDmin = 0.109 * Minidie + 0.394;
                Minimax = (Minimax > MWDmin) ? (Minimax - 0.394)/0.109 : Minimax;
                Minimin = (Minimin > MMSmax) ? (Minimin - 0.394)/0.109
                        : (BINARY == 1) ? 0.5 * Minimin : Minimin;
            } else {
                Minimin = (BINARY == 1) ? 0.5 * Minimin : Minimin;
            }
            nbinMmin = (int)floor((log10(Minimin) - ctx.stellar.logMst) / ctx.stellar.dlogM);
            nbinMmax = (int)floor((log10(Minimax) - ctx.stellar.logMst) / ctx.stellar.dlogM) + 1;
            if (nbinMmin < 0)   nbinMmin = 0;
            if (nbinMmax > ctx.stellar.nm)  nbinMmax = ctx.stellar.nm;
        }

        double logM, M_l;
        if (nbinMmax - nbinMmin > 0) {
            ran = ctx.runtime.rng->uniform() * (PlogM_cum_norm_B[nbinMmax] - PlogM_cum_norm_B[nbinMmin])
                 + PlogM_cum_norm_B[nbinMmin];
            addGammaIS *= PlogM_cum_norm_B[nbinMmax] - PlogM_cum_norm_B[nbinMmin];
        } else {
            ran = ctx.runtime.rng->uniform();
        }
        inttmp = (int)(ran * 20);
        kst = 1;
        for (int itmp = inttmp; itmp > 0; itmp--) {
            kst = imptiles_B[itmp] - 1;
            if (kst > 0) break;
        }
        logM = getcumu2xist(ctx.stellar.nm, logMass_B, PlogM_cum_norm_B, PlogM_B, ran, kst, 0);
        M_l  = pow(10.0, logM);
        double Morg = M_l;

        // --- Remnant evolution ---
        int fREM = 0;
        if (M_l > Minidie) {
            if (REMNANT == 1 || onlyWD == 1) {
                double pout[2] = {};
                Mini2Mrem(ctx, pout, M_l, 0);
                M_l  = pout[0];
                fREM = (int)pout[1];
                if (fREM >= 2) {
                    if (onlyWD == 1) { j--; continue; }
                    murell = murells[fREM - 1];
                    murelb = murelbs[fREM - 1];
                    murel  = murels[fREM - 1];
                }
            } else { j--; continue; }
        }

        // --- Binary system ---
        int swl = 0;
        double M_l2 = 99999, q21 = 99999, u0S = 99999;
        double al = 99999, alpmin = 99999, apdetS = 99999, apdetL = 99999;
        if (BINARY && fREM == 0) {
            double mult = 0.196 + 0.255*M_l;
            if (mult > MAXMULT) mult = MAXMULT;
            ran = ctx.runtime.rng->uniform();
            double coeff, q2;
            if (ran < 0.5 * mult) {
                swl = 1;
                double gamma = 1.16 - 2.79*log10(M_l);
                if (gamma > MAXGAMMA) gamma = MAXGAMMA;
                if (gamma < MINGAMMA) gamma = MINGAMMA;
                coeff = -1;
                double tmp = pow(0.1, gamma+1);
                q2 = pow((1-tmp)*ctx.runtime.rng->uniform() + tmp, 1.0/(gamma+1));
            } else if (ran < mult) {
                swl = 2;
                double gamma = (M_l >= 0.344) ? 0 : -3.09 - 6.67*log10(M_l);
                if (gamma > MAXGAMMA) gamma = MAXGAMMA;
                if (gamma < MINGAMMA) gamma = MINGAMMA;
                coeff = 1;
                double tmp = pow(0.1, gamma+1);
                q2 = pow((1-tmp)*ctx.runtime.rng->uniform() + tmp, 1.0/(gamma+1));
            } else { coeff = 0; q2 = 0; }
            if (swl > 0) {
                M_l2 = M_l * q2;
                double Ptmp = sqrt(M_l) / (sqrt(M_l) + sqrt(M_l2));
                ran = ctx.runtime.rng->uniform();
                double Mtmp1 = (ran < Ptmp) ? M_l  : M_l2;
                double Mtmp2 = (ran < Ptmp) ? M_l2 : M_l;
                double qtmp  = Mtmp2 / Mtmp1;
                double thetaE1  = sqrt(Mtmp1 * pirel * KAPPA);
                double thetaE12 = sqrt((Mtmp1+Mtmp2) * pirel * KAPPA);
                double pout[2] = {};
                getaproj(ctx, pout, M_l, M_l2, (int)coeff);
                al     = (pout[0] < 99) ? pow(10.0, pout[0]) : -1;
                alpmin = pout[1];
                u0S    = (u0obs > 0) ? u0obs : ctx.runtime.rng->uniform();
                double tmp = qtmp / u0S;
                apdetL = (sqrt(tmp+1) + sqrt(tmp)) * thetaE1 * D_l * 0.001;
                if (qtmp > 1) tmp = 1.0/qtmp/u0S;
                apdetS = (sqrt(tmp+1) - sqrt(tmp)) * thetaE12 * D_l * 0.001;
                apdetS *= (qtmp > 1) ? sqrt(1+1.0/qtmp) : sqrt(1+qtmp);
                if (alpmin > apdetL) {
                    swl += 10; M_l = Mtmp1; M_l2 = Mtmp2; q21 = qtmp;
                } else if (alpmin < apdetS) {
                    swl += 20; M_l = Mtmp1 + Mtmp2; M_l2 = 99; q21 = qtmp;
                } else { j--; continue; }
            }
        }

        // --- Microlensing observables ---
        double thetaE = sqrt(M_l * pirel * KAPPA);
        if (thetaE == 0) { printf("#ERROR: thetaE = 0!! M_l=  %.5e\n", M_l); j--; continue; }
        double tE   = thetaE / murel * 365.25;
        double piE  = pirel / thetaE;
        double piEN = piE * ( murelb*cosPA + murell*sinPA) / murel;
        double piEE = piE * (-murelb*sinPA + murell*cosPA) / murel;
        double Gamma = 1.6e-08 * D_l*D_l * thetaE * murel;
        if (NoGAMMAIS == 1) {
            Gamma = 8e-09 * D_l*D_l * thetaE * murel;
            Gamma *= addGammaIS;
            addGammaIS = 1;
        }

        // BH kick Gamma correction
        int ien = (BHhb == 1) ? 9 : 8;
        if (BHhd == 1 && i_l < ien && fREM == 3)
            Gamma *= grid.fBH_at_distance(i_l, D_l);

        stats.SumGamma += Gamma * addGammaIS;
        stats.SumtE    += Gamma * addGammaIS * tE;
        int ilogtE = (int)((log10(tE) - stats.logtEmin) /
                           (stats.logtEmax - stats.logtEmin) * stats.NbintE);
        if (ilogtE < 0) ilogtE = 0;
        if (ilogtE > stats.NbintE - 1) ilogtE = stats.NbintE - 1;
        stats.NlogtEs[ilogtE] += Gamma * addGammaIS;

        if (Gamma < ctx.runtime.rng->uniform() && SMALLGAMMA == 0) { j--; continue; }

        double addGamma = 1.0;
        double wtj = (Gamma > 1 || SMALLGAMMA == 1) ? Gamma : 1.0;

        if (wtM_L != 0) addGamma /= pow(Morg / 0.1, wtM_L);
        if (wtD_L != 0) addGamma /= pow((D_l + 1000) / 4500.0, wtD_L);

        // --- Observational likelihood ---
        double Gamma_tE     = like_obs(ctx, tE, tEobs, tEe, fetE, tEdet, UNIFORM);
        int    like_tE      = (Gamma_tE > 0) ? 1 : 0;
        if (Gamma_tE > 0) addGamma *= Gamma_tE;

        double Gamma_thetaE = like_obs(ctx, thetaE, thetaEobs, thetaEe, fethetaE, thetaEdet, UNIFORM);
        int    like_thetaE  = (Gamma_thetaE > 0) ? 1 : 0;
        if (Gamma_thetaE > 0) addGamma *= Gamma_thetaE;

        int like_piE = 1;
        if (piEe > 0) {
            double Gamma_piE = like_obs(ctx, piE, piEobs, piEe, fepiE, piEdet, UNIFORM);
            like_piE = (Gamma_piE > 0) ? 1 : 0;
            if (Gamma_piE > 0) addGamma *= Gamma_piE;
        } else if (piENe > 0 && piEEe > 0) {
            double Gamma_piEN = like_obs(ctx, piEN, piENobs, piENe, fepiEN, 0, UNIFORM);
            int like_piEN = (Gamma_piEN > 0) ? 1 : 0;
            if (Gamma_piEN > 0) addGamma *= Gamma_piEN;
            double Gamma_piEE = like_obs(ctx, piEE, piEEobs, piEEe, fepiEE, 0, UNIFORM);
            int like_piEE = (Gamma_piEE > 0) ? 1 : 0;
            if (Gamma_piEE > 0) addGamma *= Gamma_piEE;
            like_piE = like_piEN * like_piEE;
        }

        int like_mus = 1;
        if (musRCG == 1) {
            double musfac = 1.0 - D_s / Dmean;
            double vxr = vx_s - musfac*ctx.kinematics.vxsun;
            double vyr = vy_s - musfac*ctx.kinematics.vysun;
            double vzr = vz_s - musfac*ctx.kinematics.vzsun;
            muSl = (vxr*sinl      + vyr*cosl)*KS2MY/D_s;
            muSb = (vxr*cosl*sinb - vyr*sinl*sinb + vzr*cosb)*KS2MY/D_s;
        }
        if (musle > 0 && musbe > 0) {
            double G_musl = like_obs(ctx, muSl, muslobs, musle, femusl, 0, UNIFORM);
            double G_musb = like_obs(ctx, muSb, musbobs, musbe, femusb, 0, UNIFORM);
            like_mus = ((G_musl > 0) ? 1 : 0) * ((G_musb > 0) ? 1 : 0);
            if (G_musl > 0) addGamma *= G_musl;
            if (G_musb > 0) addGamma *= G_musb;
        } else if (musNe > 0 && musEe > 0) {
            double muSN =  muSb*cosPA + muSl*sinPA;
            double muSE = -muSb*sinPA + muSl*cosPA;
            double G_musN = like_obs(ctx, muSN, musNobs, musNe, femusN, 0, UNIFORM);
            double G_musE = like_obs(ctx, muSE, musEobs, musEe, femusE, 0, UNIFORM);
            like_mus = ((G_musN > 0) ? 1 : 0) * ((G_musE > 0) ? 1 : 0);
            if (G_musN > 0) addGamma *= G_musN;
            if (G_musE > 0) addGamma *= G_musE;
        }

        int like_muhel = 1;
        double muhelN = 0, muhelE = 0;
        if ((muhelNe > 0 && muhelEe > 0) || VERBOSITY >= 4) {
            double murellhel = murell + vEarthl*KS2MY*(D_s-D_l)/D_s/D_l;
            double murelbhel = murelb + vEarthb*KS2MY*(D_s-D_l)/D_s/D_l;
            muhelN =  murelbhel*cosPA + murellhel*sinPA;
            muhelE = -murelbhel*sinPA + murellhel*cosPA;
        }
        if (muhelNe > 0 && muhelEe > 0) {
            double G_N = like_obs(ctx, muhelN, muhelNobs, muhelNe, femuhelN, 0, UNIFORM);
            double G_E = like_obs(ctx, muhelE, muhelEobs, muhelEe, femuhelE, 0, UNIFORM);
            like_muhel = ((G_N > 0) ? 1 : 0) * ((G_E > 0) ? 1 : 0);
            if (G_N > 0) addGamma *= G_N;
            if (G_E > 0) addGamma *= G_E;
        }

        int like_IL = 1;
        double IL = 99;
        if (AI0 > 0 && fREM == 0) {
            double M1 = (swl < 20) ? M_l : M_l / (1 + q21);
            IL = mf.I_magnitude(M1);
            if (swl >= 20) {
                double M2  = q21 * M_l / (1 + q21);
                double IL2 = mf.I_magnitude(M2);
                double F12 = pow(10, -0.4*IL) + pow(10, -0.4*IL2);
                IL = -2.5 * log10(F12);
            }
            IL += ext->distance_modulus_term(D_l) + ext->at_distance(D_l).i_band;
            if (ILe > 0) {
                double G_IL = like_obs(ctx, IL, ILobs, ILe, feIL, ILdet, UNIFORM);
                like_IL = (G_IL > 0) ? 1 : 0;
                if (G_IL > 0) addGamma *= G_IL;
            }
        }

        int like_KL = 1;
        double KL = 99;
        if (AK0 > 0 && fREM == 0) {
            double M1 = (swl < 20) ? M_l : M_l / (1 + q21);
            KL = mf.K_magnitude(M1);
            if (swl >= 20) {
                double M2  = q21 * M_l / (1 + q21);
                double KL2 = mf.K_magnitude(M2);
                double F12 = pow(10, -0.4*KL) + pow(10, -0.4*KL2);
                KL = -2.5 * log10(F12);
            }
            KL += ext->distance_modulus_term(D_l) + ext->at_distance(D_l).k_band;
            if (KLe > 0) {
                double G_KL = like_obs(ctx, KL, KLobs, KLe, feKL, KLdet, UNIFORM);
                like_KL = (G_KL > 0) ? 1 : 0;
                if (G_KL > 0) addGamma *= G_KL;
            }
        }

	        // Optional custom likelihood
	        Event event;
	        event.weight = wtj;
	        event.tE = tE;
	        event.thetaE = thetaE;
	        event.piE = piE;
	        event.piEN = piEN;
	        event.piEE = piEE;
	        event.lens_distance_pc = D_l;
	        event.source_distance_pc = D_s;
	        event.lens_mass_msun = M_l;
	        event.mu_rel_masyr = murel;
            event.mu_rel_N_masyr = murelb*cosPA + murell*sinPA;
            event.mu_rel_E_masyr = -murelb*sinPA + murell*cosPA;
            event.source_mu_l_masyr = muSl;
            event.source_mu_b_masyr = muSb;
            event.lens_i_mag = IL;
            event.lens_k_mag = KL;
	        event.lens_component = i_l;
	        event.source_component = i_s;
            event.remnant_flag = fREM;
            event.source_age_gyr = tau_s;
            event.lens_age_gyr = tau_l;
            if (attach_source_properties) {
                try {
                    model::ForwardSourceQuery source_query;
                    // The current PARSEC source tables do not include a
                    // stellar-halo sequence. Keep the event component label as
                    // halo, but use the old, metal-poor thick-disk sequence as
                    // a conservative source-property proxy.
                    source_query.component_index = (i_s == 10) ? 7 : i_s;
                    source_query.distance_pc = D_s;
                    source_query.min_initial_mass_msun = cfg.source_min_initial_mass_msun;
                    source_query.max_initial_mass_msun = cfg.source_max_initial_mass_msun;
                    const auto population_point = sample_source_population_point(
                        *cfg.forward_source_generator, i_s, D_s, cfg);
                    source_query.use_default_log_age = false;
                    source_query.log_age = population_point.log_age;
                    source_query.use_default_metallicity = false;
                    source_query.metallicity_mh = population_point.metallicity_mh;
                    const auto source = sample_matched_source(*cfg.forward_source_generator,
                                                              source_query,
                                                              cfg);
                    event.source_log_age = source.stellar.log_age;
                    event.source_metallicity_mh = source.stellar.metallicity_mh;
                    event.source_zini = source.stellar.zini;
                    event.source_initial_mass_msun = source.stellar.initial_mass_msun;
                    event.source_current_mass_msun = source.stellar.current_mass_msun;
                    event.source_radius_rsun = source.stellar.radius_rsun;
                    event.source_teff_k = source.stellar.teff_k;
                    event.source_logg = source.stellar.logg;
                    event.source_angular_radius_microarcsec = source.angular_radius_microarcsec;
                    for (const auto &band : cfg.forward_source_generator->bands()) {
                        const auto found = source.stellar.absolute_magnitudes.find(band);
                        event.source_absolute_magnitudes.push_back(
                            found == source.stellar.absolute_magnitudes.end()
                                ? std::numeric_limits<double>::quiet_NaN()
                                : found->second);
                    }
                } catch (const std::exception &) {
                    j--;
                    continue;
                }
            }
	        if (custom_likelihood) {
	            double cl = custom_likelihood(event);
	            if (cl <= 0) { j--; continue; }
	            addGamma *= cl;
	        }

        Gamma *= addGamma * addGammaIS;
        wtj   *= addGamma * addGammaIS;

        int like     = like_tE * like_thetaE * like_piE * like_mus * like_muhel * like_IL * like_KL;
        int like_expiE = like_tE * like_thetaE * like_mus * like_muhel * like_IL * like_KL;
        int like_exthE = like_tE * like_piE    * like_mus * like_muhel * like_IL * like_KL;

        if (CALCPRIORpiE && piENe > 0 && piEEe > 0 && like_expiE == 1) {
            double chi2N = (piEN - piENobs)/piENe; chi2N *= chi2N;
            double chi2E = (piEE - piEEobs)/piEEe; chi2E *= chi2E;
            stats.wtlike_except_piE += wtj;
            stats.wtlike_w_piEe     += wtj * 0.5/PI/piENe/piEEe * exp(-0.5*(chi2N+chi2E));
        }
        if (CALCPRIORthE && thetaEe > 0 && like_exthE == 1) {
            double chi2t = (thetaE - thetaEobs)/thetaEe; chi2t *= chi2t;
            stats.wtlike_except_thE += wtj;
            stats.wtlike_w_thEe     += wtj * sqrt(0.5/PI)/thetaEe * exp(-0.5*chi2t);
        }

	        if (like_tE == 1) stats.wtlike_tE += wtj;
	        if (like == 1) { stats.Nlike++; stats.wtlike += wtj; }
	        else continue;
	        event.weight = wtj;
	        if (event_sink) event_sink(event);

	        // --- Output ---
	        EventOutputRow output;
	        output.wtj = wtj;
	        output.tE = tE;
	        output.thetaE = thetaE;
	        output.piE = piE;
	        output.piEN = piEN;
	        output.piEE = piEE;
	        output.lens_mass = M_l;
	        output.lens_distance = D_l;
	        output.source_distance = D_s;
	        output.vx_s = vx_s;
	        output.vy_s = vy_s;
	        output.vz_s = vz_s;
	        output.vx_l = vx_l;
	        output.vy_l = vy_l;
	        output.vz_l = vz_l;
	        output.mu_rel = murel;
	        output.muSl = muSl;
	        output.muSb = muSb;
	        output.IL = IL;
	        output.KL = KL;
	        output.source_component = i_s;
	        output.lens_component = i_l;
	        output.tau_s = tau_s;
	        output.tau_l = tau_l;
	        output.remnant_flag = fREM;
	        output.gamma = Gamma;
	        output.muhelN = muhelN;
	        output.muhelE = muhelE;
	        output.thetaE_lens_radius = thetaE * D_l * 1e-03;
	        output.q21 = q21;
	        output.secondary_mass = M_l2;
	        output.projected_separation = al;
	        output.projected_separation_min = alpmin;
	        output.u0_source = u0S;
	        output.binary_state = swl;
	        reporter.print_event(output);

        // Count remnants and binary types
        stats.ncntall += wtj;
        if (swl == 0)               stats.ncnts   += wtj;
        if (swl > 10 && swl < 20)   stats.ncntbWD += wtj;
        if (swl > 20)               stats.ncntbCD += wtj;
	        if (emit_cli_output && apdetS > apdetL)
	            printf("# ERROR: apdetS (= %.6f) > apdetL (= %.6f) !!!!!\n", apdetS, apdetL);
        if (fREM == 0 && M_l < 0.08) stats.nBD += wtj;
        if (fREM == 0 && M_l > 0.08) stats.nMS += wtj;
        if (fREM == 1) stats.nWD += wtj;
        if (fREM == 2) stats.nNS += wtj;
        if (fREM == 3) stats.nBH += wtj;
    }

    // --- Summary ---
    if (cfg.rate_summary != nullptr) {
        *cfg.rate_summary = make_rate_summary(stats, NSIMU, cfg.l, cfg.b,
                                              cfg.Nsall, cfg.nallS, cfg.tauall);
    }
    reporter.print_summary(stats, NSIMU, CALCPRIORpiE, CALCPRIORthE,
                           cfg.Nsall, cfg.nallS, cfg.tauall);

	    return 0;
	}

int EventSampler::run_cli(RunContext &ctx,
                          LineOfSightDensityGrid &grid,
                          PopulationRuntime &pop,
                          MassFunction &mf,
                          const Config &cfg,
                          const ObservationConfig &obs)
{
    return run(ctx, grid, pop, mf, cfg, obs, {}, {}, true);
}

} // namespace genulens
