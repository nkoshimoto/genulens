#pragma once

#include <filesystem>
#include <cstdio>
#include <string>
#include <vector>

namespace genulens {

class InputDataRepository {
public:
    explicit InputDataRepository(std::filesystem::path override_dir = {});

    std::filesystem::path resolve(const std::string &filename) const;
    std::vector<std::filesystem::path> search_roots() const;

private:
    std::filesystem::path override_dir_;
};

std::filesystem::path resolve_input_file(const std::string &filename);
FILE *open_input_file(const char *filename, const char *mode);
const char *resolve_input_path_c(const char *filename);

} // namespace genulens
