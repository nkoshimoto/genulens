#include "genulens/options.hpp"

#include <cstdlib>
#include <cstring>

namespace genulens {

namespace {

const char *get_option(int argc, char **argv, const char *name, int argno)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], name) == 0 && i + argno < argc) {
            return argv[i + argno];
        }
    }
    return nullptr;
}

} // namespace

GenulensConfig config_from_cli(int argc, char **argv)
{
    GenulensConfig cfg;
    cfg.raw_cli_args.assign(argv, argv + argc);
    if (const char *value = get_option(argc, argv, "seed", 1)) cfg.seed = std::strtoul(value, nullptr, 10);
    if (const char *value = get_option(argc, argv, "l", 1)) cfg.l = std::atof(value);
    if (const char *value = get_option(argc, argv, "b", 1)) cfg.b = std::atof(value);
    if (const char *value = get_option(argc, argv, "Nsimu", 1)) cfg.n_simu = std::atol(value);
    if (const char *value = get_option(argc, argv, "tE", 1)) cfg.observed_tE = std::atof(value);
    if (const char *value = get_option(argc, argv, "tE", 2)) cfg.observed_tE_error = std::atof(value);
    return cfg;
}

} // namespace genulens

