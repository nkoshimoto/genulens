#pragma once

#include <filesystem>
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

} // namespace genulens

