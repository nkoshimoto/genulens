#include "genulens/model/isochrone_library.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace genulens::model {
namespace {

std::string alpha_token(double alpha_fe)
{
    const char prefix = alpha_fe >= 0.0 ? 'p' : 'm';
    const double rounded = std::round(std::abs(alpha_fe) * 100.0) / 100.0;
    std::ostringstream out;
    out.setf(std::ios::fixed);
    out.precision(2);
    out << "afe_" << prefix << rounded;
    auto token = out.str();
    for (auto &ch : token) {
        if (ch == '.') ch = 'p';
    }
    return token;
}

void validate_photometry(const std::string &photometry)
{
    if (photometry != "roman" && photometry != "prime") {
        throw std::runtime_error("unknown forward source photometry: " + photometry);
    }
}

} // namespace

std::string normalize_isochrone_family(const std::string &family)
{
    if (family.empty() || family == "parsec" || family == "parsec_cmd") return "parsec";
    if (family == "mist" || family == "MIST") return "mist";
    throw std::runtime_error(
        "unknown isochrone family: " + family + " (expected parsec or mist)");
}

std::string normalize_isochrone_abundance(const std::string &abundance, double alpha_fe)
{
    if (abundance.empty() || abundance == "solar" || abundance == "solar_scaled") {
        return "solar_scaled";
    }
    if (abundance == "alpha" || abundance == "alpha_enhanced") {
        return "alpha_enhanced";
    }
    if (std::abs(alpha_fe) > 0.0) return "alpha_enhanced";
    throw std::runtime_error(
        "unknown isochrone abundance: " + abundance +
        " (expected solar_scaled or alpha_enhanced)");
}

std::string default_isochrone_table_path(const IsochroneLibrarySpec &input)
{
    const auto family = normalize_isochrone_family(input.family);
    const auto abundance = normalize_isochrone_abundance(input.abundance, input.alpha_fe);
    validate_photometry(input.photometry);

    if (family == "parsec" && abundance == "solar_scaled") {
        return "source_photometry/parsec_cmd/metallicity_grid/normalized/all_" +
               input.photometry + "_parsec.dat";
    }

    if (family == "parsec" && abundance == "alpha_enhanced") {
        throw std::runtime_error(
            "PARSEC/CMD alpha-enhanced source-photometry tables are not available in "
            "genulens. Use isochrone_family='mist' for alpha-enhanced "
            "experiments, or set isochrone_abundance='solar_scaled' for PARSEC.");
    }

    if (family == "mist" && abundance == "solar_scaled") {
        return "source_photometry/mist/v2.5/solar_scaled/normalized/all_" +
               input.photometry + "_mist_solar.dat";
    }

    if (family == "mist" && abundance == "alpha_enhanced") {
        return "source_photometry/mist/v2.5/alpha_enhanced/" + alpha_token(input.alpha_fe) +
               "/normalized/all_" + input.photometry + "_mist_alpha_" +
               alpha_token(input.alpha_fe) + ".dat";
    }

    throw std::runtime_error("unsupported isochrone library specification");
}

} // namespace genulens::model
