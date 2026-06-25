#include "genulens/io/input_data.hpp"

#include <cstdlib>
#include <stdexcept>
#include <string>

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
#ifdef GENULENS_SOURCE_DIR
    roots.emplace_back(std::filesystem::path(GENULENS_SOURCE_DIR) / "input_files");
#endif
    roots.emplace_back("/usr/local/share/genulens/input_files");
    roots.emplace_back("/usr/share/genulens/input_files");
    return roots;
}

std::filesystem::path InputDataRepository::resolve(const std::string &filename) const
{
    const auto direct = std::filesystem::path(filename);
    if (std::filesystem::exists(direct)) {
        return direct;
    }

    auto relative = direct;
    if (!relative.empty() && *relative.begin() == "input_files") {
        std::filesystem::path stripped;
        auto it = relative.begin();
        ++it;
        for (; it != relative.end(); ++it) stripped /= *it;
        relative = stripped;
    }

    for (const auto &root : search_roots()) {
        const auto candidate = root / relative;
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    const auto basename = direct.filename();
    for (const auto &root : search_roots()) {
        const auto candidate = root / basename;
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    throw std::runtime_error("could not resolve input file: " + filename);
}

std::filesystem::path resolve_input_file(const std::string &filename)
{
    return InputDataRepository().resolve(filename);
}

FILE *open_input_file(const char *filename, const char *mode)
{
    if (!filename) return nullptr;
    if (FILE *fp = std::fopen(filename, mode)) {
        return fp;
    }
    try {
        const auto resolved = resolve_input_file(filename);
        return std::fopen(resolved.string().c_str(), mode);
    } catch (const std::exception &) {
        return nullptr;
    }
}

const char *resolve_input_path_c(const char *filename)
{
    static thread_local std::string resolved;
    resolved = resolve_input_file(filename).string();
    return resolved.c_str();
}

} // namespace genulens
