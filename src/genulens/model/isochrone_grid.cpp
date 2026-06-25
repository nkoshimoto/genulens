#include "genulens/model/isochrone_grid.hpp"

#include "genulens/io/input_data.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace genulens::model {
namespace {

std::vector<std::string> split_ws(const std::string &line)
{
    std::istringstream in(line);
    std::vector<std::string> out;
    std::string token;
    while (in >> token) out.push_back(token);
    return out;
}

double parse_double(const std::string &value, const std::filesystem::path &path)
{
    try {
        return std::stod(value);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid numeric value in " + path.string() + ": " + value);
    }
}

int parse_int(const std::string &value, const std::filesystem::path &path)
{
    try {
        return std::stoi(value);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid integer value in " + path.string() + ": " + value);
    }
}

double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}

template <typename Row>
bool continuous_isochrone_segment(const Row &lo, const Row &hi)
{
    if (!(hi.initial_mass_msun > lo.initial_mass_msun)) return false;
    if (!(lo.teff_k > 0.0 && hi.teff_k > 0.0)) return false;
    if (!(lo.radius_rsun > 0.0 && hi.radius_rsun > 0.0)) return false;

    const double dlog_teff = std::abs(std::log10(hi.teff_k) - std::log10(lo.teff_k));
    const double dlog_radius = std::abs(std::log10(hi.radius_rsun) - std::log10(lo.radius_rsun));
    const double dlogg = std::abs(hi.logg - lo.logg);
    if (dlog_teff > 0.12 || dlog_radius > 0.45 || dlogg > 0.75) return false;

    const std::size_t nmag = std::min(lo.magnitudes.size(), hi.magnitudes.size());
    for (std::size_t i = 0; i < nmag; ++i) {
        if (!std::isfinite(lo.magnitudes[i]) || !std::isfinite(hi.magnitudes[i])) {
            return false;
        }
        if (std::abs(hi.magnitudes[i] - lo.magnitudes[i]) > 2.0) {
            return false;
        }
    }
    return true;
}

void intersect_linear_constraint(double y0, double y1, double ymin, double ymax,
                                 double &tlo, double &thi)
{
    const double dy = y1 - y0;
    if (dy == 0.0) {
        if (!(y0 >= ymin && y0 <= ymax)) {
            tlo = 1.0;
            thi = 0.0;
        }
        return;
    }

    double a = (ymin - y0) / dy;
    double b = (ymax - y0) / dy;
    if (a > b) std::swap(a, b);
    tlo = std::max(tlo, a);
    thi = std::min(thi, b);
}

} // namespace

IsochroneGrid IsochroneGrid::load_default_roman()
{
    return load("source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat");
}

IsochroneGrid IsochroneGrid::load_default_prime()
{
    return load("source_photometry/parsec_cmd/metallicity_grid/normalized/all_prime_parsec.dat");
}

IsochroneGrid IsochroneGrid::load(const std::filesystem::path &path)
{
    const auto resolved = genulens::resolve_input_file(path.string());
    std::ifstream in(resolved);
    if (!in) {
        throw std::runtime_error("could not open isochrone grid: " + resolved.string());
    }

    IsochroneGrid grid;
    std::vector<std::string> columns;
    std::unordered_map<std::string, std::size_t> index;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            columns = split_ws(line.substr(1));
            index.clear();
            for (std::size_t i = 0; i < columns.size(); ++i) index[columns[i]] = i;
            grid.bands_.clear();
            const auto logg_column = index.find("logg");
            if (logg_column == index.end()) {
                throw std::runtime_error("missing isochrone column: logg");
            }
            for (std::size_t i = logg_column->second + 1; i < columns.size(); ++i) {
                grid.bands_.push_back(columns[i]);
            }
            continue;
        }

        if (columns.empty()) {
            throw std::runtime_error("isochrone grid missing header: " + resolved.string());
        }
        const auto values = split_ws(line);
        if (values.size() != columns.size()) {
            throw std::runtime_error("isochrone grid row has wrong column count in " + resolved.string());
        }

        auto at = [&](const std::string &name) -> const std::string & {
            const auto found = index.find(name);
            if (found == index.end()) throw std::runtime_error("missing isochrone column: " + name);
            return values[found->second];
        };

        Row row;
        row.component = at("component");
        row.component_index = parse_int(at("component_index"), resolved);
        row.family = at("family");
        row.log_age = parse_double(at("logAge"), resolved);
        row.metallicity_mh = parse_double(at("MH"), resolved);
        row.zini = parse_double(at("Zini"), resolved);
        row.initial_mass_msun = parse_double(at("Mini"), resolved);
        row.current_mass_msun = parse_double(at("Mass"), resolved);
        row.radius_rsun = parse_double(at("radius_rsun"), resolved);
        row.teff_k = parse_double(at("Teff_K"), resolved);
        row.logg = parse_double(at("logg"), resolved);
        row.magnitudes.reserve(grid.bands_.size());
        for (const auto &band : grid.bands_) {
            row.magnitudes.push_back(parse_double(at(band), resolved));
        }
        grid.rows_.push_back(std::move(row));
    }

    if (grid.rows_.empty()) {
        throw std::runtime_error("isochrone grid has no rows: " + resolved.string());
    }
    if (grid.bands_.empty()) {
        throw std::runtime_error("isochrone grid has no magnitude columns: " + resolved.string());
    }
    grid.build_sequences();
    return grid;
}

void IsochroneGrid::build_sequences()
{
    std::stable_sort(rows_.begin(), rows_.end(), [](const Row &a, const Row &b) {
        if (a.component != b.component) return a.component < b.component;
        if (a.log_age != b.log_age) return a.log_age < b.log_age;
        if (a.metallicity_mh != b.metallicity_mh) return a.metallicity_mh < b.metallicity_mh;
        return a.initial_mass_msun < b.initial_mass_msun;
    });

    sequences_.clear();
    for (const auto &row : rows_) {
        if (sequences_.empty() ||
            sequences_.back().component != row.component ||
            sequences_.back().log_age != row.log_age ||
            sequences_.back().metallicity_mh != row.metallicity_mh) {
            Sequence seq;
            seq.component = row.component;
            seq.family = row.family;
            seq.component_index = row.component_index;
            seq.log_age = row.log_age;
            seq.metallicity_mh = row.metallicity_mh;
            sequences_.push_back(std::move(seq));
        }
        auto &rows = sequences_.back().rows;
        if (!rows.empty() &&
            std::abs(rows.back().initial_mass_msun - row.initial_mass_msun) < 1.0e-12) {
            continue;
        }
        rows.push_back(row);
    }
}

const IsochroneGrid::Sequence &IsochroneGrid::select_sequence(const IsochroneQuery &query) const
{
    const Sequence *best = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const auto &seq : sequences_) {
        if (!query.component.empty() && seq.component != query.component) continue;
        const double age_score = std::abs(seq.log_age - query.log_age);
        const double mh_score = std::abs(seq.metallicity_mh - query.metallicity_mh);
        const double score = age_score * 100.0 + mh_score;
        if (score < best_score) {
            best = &seq;
            best_score = score;
        }
    }
    if (!best) {
        throw std::runtime_error("no isochrone sequence matches query component: " + query.component);
    }
    return *best;
}

StellarProperties IsochroneGrid::lookup(const IsochroneQuery &query) const
{
    const auto &seq = select_sequence(query);
    if (seq.rows.empty()) throw std::runtime_error("selected empty isochrone sequence");

    const Row *lo = &seq.rows.front();
    const Row *hi = &seq.rows.back();
    if (query.initial_mass_msun <= seq.rows.front().initial_mass_msun) {
        hi = lo;
    } else if (query.initial_mass_msun >= seq.rows.back().initial_mass_msun) {
        lo = hi;
    } else {
        auto upper = std::lower_bound(seq.rows.begin(), seq.rows.end(), query.initial_mass_msun,
                                      [](const Row &row, double mass) {
                                          return row.initial_mass_msun < mass;
                                      });
        hi = &(*upper);
        lo = &(*(upper - 1));
    }

    double t = 0.0;
    if (hi->initial_mass_msun != lo->initial_mass_msun) {
        t = (query.initial_mass_msun - lo->initial_mass_msun) /
            (hi->initial_mass_msun - lo->initial_mass_msun);
    }
    if (!continuous_isochrone_segment(*lo, *hi)) {
        t = (t < 0.5) ? 0.0 : 1.0;
    }

    StellarProperties props;
    props.component = lo->component;
    props.component_index = lo->component_index;
    props.family = lo->family;
    props.log_age = lerp(lo->log_age, hi->log_age, t);
    props.metallicity_mh = lerp(lo->metallicity_mh, hi->metallicity_mh, t);
    props.zini = lerp(lo->zini, hi->zini, t);
    props.initial_mass_msun = lerp(lo->initial_mass_msun, hi->initial_mass_msun, t);
    props.current_mass_msun = lerp(lo->current_mass_msun, hi->current_mass_msun, t);
    props.radius_rsun = lerp(lo->radius_rsun, hi->radius_rsun, t);
    props.teff_k = lerp(lo->teff_k, hi->teff_k, t);
    props.logg = lerp(lo->logg, hi->logg, t);
    for (std::size_t i = 0; i < bands_.size(); ++i) {
        props.absolute_magnitudes[bands_[i]] = lerp(lo->magnitudes[i], hi->magnitudes[i], t);
    }
    return props;
}

std::vector<MassInterval> IsochroneGrid::matching_initial_mass_intervals(
    const IsochroneQuery &query,
    const std::vector<MagnitudeSelection> &selection) const
{
    const auto &seq = select_sequence(query);
    if (seq.rows.size() < 2 || selection.empty()) return {};

    std::vector<std::size_t> band_indices;
    band_indices.reserve(selection.size());
    for (const auto &cut : selection) {
        const auto found = std::find(bands_.begin(), bands_.end(), cut.band);
        if (found == bands_.end()) {
            throw std::runtime_error("isochrone band is not available: " + cut.band);
        }
        band_indices.push_back(static_cast<std::size_t>(found - bands_.begin()));
    }

    std::vector<MassInterval> intervals;
    for (std::size_t i = 0; i + 1 < seq.rows.size(); ++i) {
        const auto &lo = seq.rows[i];
        const auto &hi = seq.rows[i + 1];
        if (!continuous_isochrone_segment(lo, hi)) continue;

        double tlo = 0.0;
        double thi = 1.0;
        for (std::size_t j = 0; j < selection.size(); ++j) {
            const auto &cut = selection[j];
            const std::size_t band_index = band_indices[j];
            const double mag0 = lo.magnitudes[band_index] + cut.magnitude_offset;
            const double mag1 = hi.magnitudes[band_index] + cut.magnitude_offset;
            if (!std::isfinite(mag0) || !std::isfinite(mag1)) {
                tlo = 1.0;
                thi = 0.0;
                break;
            }
            intersect_linear_constraint(mag0, mag1, cut.min_magnitude,
                                        cut.max_magnitude, tlo, thi);
            if (tlo > thi) break;
        }
        if (tlo > thi) continue;

        MassInterval interval;
        interval.min_mass_msun = lerp(lo.initial_mass_msun, hi.initial_mass_msun, tlo);
        interval.max_mass_msun = lerp(lo.initial_mass_msun, hi.initial_mass_msun, thi);
        if (!(interval.max_mass_msun > interval.min_mass_msun)) continue;
        if (!intervals.empty() &&
            interval.min_mass_msun <= intervals.back().max_mass_msun * (1.0 + 1e-12)) {
            intervals.back().max_mass_msun =
                std::max(intervals.back().max_mass_msun, interval.max_mass_msun);
        } else {
            intervals.push_back(interval);
        }
    }
    return intervals;
}

} // namespace genulens::model
