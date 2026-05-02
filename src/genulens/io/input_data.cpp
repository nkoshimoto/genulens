#include "genulens/io/input_data.hpp"

#include <cstdlib>
#include <stdexcept>

namespace genulens {

InputDataRepository::InputDataRepository(std::filesystem::path override_dir)
    : override_dir_(std::move(override_dir))
{
}

std::vector<std::filesystem::path> InputDataRepository::search_roots() const
{
    std::vector<std::filesystem::path> roots;
    if (!override_dir_.empty()) roots.push_back(override_dir_);
    roots.push_back(std::filesystem::current_path() / "input_files");
    if (const char *env = std::getenv("GENULENS_INPUT_DIR")) {
        roots.emplace_back(env);
    }
    roots.emplace_back("/usr/local/share/genulens/input_files");
    roots.emplace_back("/usr/share/genulens/input_files");
    return roots;
}

std::filesystem::path InputDataRepository::resolve(const std::string &filename) const
{
    for (const auto &root : search_roots()) {
        const auto candidate = root / filename;
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    throw std::runtime_error("could not resolve input file: " + filename);
}

} // namespace genulens

