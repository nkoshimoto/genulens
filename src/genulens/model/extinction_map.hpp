#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace genulens::model {

struct ExtinctionMapSample {
    double l_deg = 0.0;
    double b_deg = 0.0;
    double ejk = 0.0;
};

class GenstarsExtinctionMap {
public:
    static GenstarsExtinctionMap load(const std::filesystem::path &path);
    static GenstarsExtinctionMap load_default(int extinction_map);

    ExtinctionMapSample sample(double l_deg, double b_deg, int extinction_map) const;
    double ejk_at(double l_deg, double b_deg, int extinction_map) const;
    const std::filesystem::path &path() const { return path_; }
    std::size_t row_count() const { return rows_.size(); }

private:
    struct Row {
        double l = 0.0;
        double b = 0.0;
        double mean_ejk = 0.0;
        std::vector<double> subgrid_ejk;
    };

    std::filesystem::path path_;
    std::vector<Row> rows_;
};

} // namespace genulens::model
