#pragma once

#include <string>

namespace genulens::model {

struct IsochroneLibrarySpec {
    std::string family = "parsec";
    std::string photometry = "roman";
    std::string abundance = "solar_scaled";
    double alpha_fe = 0.0;
};

std::string normalize_isochrone_family(const std::string &family);
std::string normalize_isochrone_abundance(const std::string &abundance, double alpha_fe);
std::string default_isochrone_table_path(const IsochroneLibrarySpec &spec);

} // namespace genulens::model
