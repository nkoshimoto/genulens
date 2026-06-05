#include "genulens/model/extinction_map.hpp"

#include "genulens/io/input_data.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace genulens::model {
namespace {

constexpr double kCellSizeDeg = 0.025;
constexpr double kHalfCellSizeDeg = 0.5 * kCellSizeDeg;
constexpr double kBoundaryTolerance = 1.0e-10;

std::vector<double> parse_numbers(const std::string &line)
{
    std::istringstream stream(line);
    std::vector<double> values;
    double value = 0.0;
    while (stream >> value) values.push_back(value);
    return values;
}

bool contains(double center, double value)
{
    return center - kHalfCellSizeDeg - kBoundaryTolerance <= value &&
           value < center + kHalfCellSizeDeg + kBoundaryTolerance;
}

int nearest_subgrid_count(std::size_t n_subgrid)
{
    const auto root = static_cast<int>(std::lround(std::sqrt(static_cast<double>(n_subgrid))));
    return root > 0 && static_cast<std::size_t>(root * root) == n_subgrid ? root : 0;
}

} // namespace

GenstarsExtinctionMap GenstarsExtinctionMap::load(const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open genstars extinction map: " + path.string());
    }

    GenstarsExtinctionMap map;
    map.path_ = path;

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto values = parse_numbers(line);
        if (values.size() < 3) continue;

        Row row;
        row.l = values[0];
        row.b = values[1];
        row.mean_ejk = values[2];
        if (values.size() > 3) {
            row.subgrid_ejk.assign(values.begin() + 3, values.end());
        }
        map.rows_.push_back(std::move(row));
    }

    if (map.rows_.empty()) {
        throw std::runtime_error("genstars extinction map has no data rows: " + path.string());
    }
    return map;
}

GenstarsExtinctionMap GenstarsExtinctionMap::load_default(int extinction_map)
{
    const std::string filename = extinction_map == 0 ? "EJK_G12_S20.dat" : "EJK_G12_S20_LR.dat";
    return load(resolve_input_file(filename));
}

ExtinctionMapSample GenstarsExtinctionMap::sample(double l_deg, double b_deg, int extinction_map) const
{
    const Row *best = nullptr;
    double best_distance2 = 0.0;
    for (const auto &row : rows_) {
        if (!contains(row.l, l_deg) || !contains(row.b, b_deg)) continue;
        best = &row;
        break;
    }
    if (best == nullptr) {
        for (const auto &row : rows_) {
            const double dl = row.l - l_deg;
            const double db = row.b - b_deg;
            const double distance2 = dl * dl + db * db;
            if (best == nullptr || distance2 < best_distance2) {
                best = &row;
                best_distance2 = distance2;
            }
        }
    }
    if (best == nullptr) {
        throw std::runtime_error("genstars extinction map is empty");
    }

    double ejk = best->mean_ejk;
    if (extinction_map < 2 && !best->subgrid_ejk.empty()) {
        const int n_side = nearest_subgrid_count(best->subgrid_ejk.size());
        if (n_side > 0) {
            const double sub_size = kCellSizeDeg / n_side;
            int il = static_cast<int>(std::floor((l_deg - (best->l - kHalfCellSizeDeg)) / sub_size));
            int ib = static_cast<int>(std::floor((b_deg - (best->b - kHalfCellSizeDeg)) / sub_size));
            il = std::clamp(il, 0, n_side - 1);
            ib = std::clamp(ib, 0, n_side - 1);
            const int index = il * n_side + ib;
            ejk = best->subgrid_ejk[static_cast<std::size_t>(index)];
        }
    }

    return {best->l, best->b, ejk};
}

double GenstarsExtinctionMap::ejk_at(double l_deg, double b_deg, int extinction_map) const
{
    return sample(l_deg, b_deg, extinction_map).ejk;
}

} // namespace genulens::model
