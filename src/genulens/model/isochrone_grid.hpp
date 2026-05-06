#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace genulens::model {

struct IsochroneQuery {
    double log_age = 0.0;
    double metallicity_mh = 0.0;
    double initial_mass_msun = 0.0;
    std::string component;
};

struct StellarProperties {
    std::string component;
    int component_index = -1;
    std::string family;
    double log_age = 0.0;
    double metallicity_mh = 0.0;
    double zini = 0.0;
    double initial_mass_msun = 0.0;
    double current_mass_msun = 0.0;
    double radius_rsun = 0.0;
    double teff_k = 0.0;
    double logg = 0.0;
    std::map<std::string, double> absolute_magnitudes;
};

class IsochroneGrid {
public:
    static IsochroneGrid load(const std::filesystem::path &path);
    static IsochroneGrid load_default_roman();
    static IsochroneGrid load_default_prime();

    const std::vector<std::string> &bands() const { return bands_; }
    std::size_t row_count() const { return rows_.size(); }
    std::size_t sequence_count() const { return sequences_.size(); }

    StellarProperties lookup(const IsochroneQuery &query) const;

private:
    struct Row {
        std::string component;
        int component_index = -1;
        std::string family;
        double log_age = 0.0;
        double metallicity_mh = 0.0;
        double zini = 0.0;
        double initial_mass_msun = 0.0;
        double current_mass_msun = 0.0;
        double radius_rsun = 0.0;
        double teff_k = 0.0;
        double logg = 0.0;
        std::vector<double> magnitudes;
    };

    struct Sequence {
        std::string component;
        std::string family;
        int component_index = -1;
        double log_age = 0.0;
        double metallicity_mh = 0.0;
        std::vector<Row> rows;
    };

    std::vector<std::string> bands_;
    std::vector<Row> rows_;
    std::vector<Sequence> sequences_;

    void build_sequences();
    const Sequence &select_sequence(const IsochroneQuery &query) const;
};

} // namespace genulens::model
