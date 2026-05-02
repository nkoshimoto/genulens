#pragma once

#include <string>
#include <vector>

namespace genulens {

struct GenulensConfig {
    unsigned long seed = 12304357UL;
    double l = 1.0;
    double b = -3.9;
    long n_simu = 100000;
    double observed_tE = 54.5;
    double observed_tE_error = 99999999999.0;
    bool use_gaussian_likelihood = true;
    std::string input_dir;
    std::vector<std::string> raw_cli_args;
};

GenulensConfig config_from_cli(int argc, char **argv);

} // namespace genulens

