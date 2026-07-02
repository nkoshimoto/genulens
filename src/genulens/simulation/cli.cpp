#include "genulens/simulation/cli.hpp"

#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/sampler.hpp"

#include <cstring>
#include <iostream>

namespace genulens {

namespace {

bool wants_help(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            return true;
        }
    }
    return false;
}

void print_help()
{
    std::cout
        << "Usage: ./genulens [name value]...\n"
        << "Example: ./genulens l 0.5 b 0.2 Nsimu 1000 seed 1234\n"
        << "Common options:\n"
        << "  l <deg>          Galactic longitude\n"
        << "  b <deg>          Galactic latitude\n"
        << "  Nsimu <count>    Number of Monte Carlo samples\n"
        << "  seed <int>       Random seed\n"
        << "  tE <days>        Observed Einstein timescale\n"
        << "  thetaE <mas>     Observed angular Einstein radius\n"
        << "  piE <value>      Observed microlens parallax\n"
        << "  IL <mag>         Lens/source brightness constraint\n"
        << "  REMNANT <0|1>    Include remnants\n"
        << "  BINARY <0|1>     Include binaries\n"
        << "Input data search order:\n"
        << "  1. direct path\n"
        << "  2. ./input_files\n"
        << "  3. GENULENS_INPUT_DIR\n"
        << "  4. installed shared data near the executable/module\n"
        << "  5. system shared data\n";
}

} // namespace

int run_cli(int argc, char **argv)
{
    if (wants_help(argc, argv)) {
        print_help();
        return 0;
    }

    Initializer initializer;
    auto context = initializer.create_context();
    initializer.initialize_rng(context, argc, argv);
    Sampler sampler;
    return sampler.run_cli(context, argc, argv);
}

} // namespace genulens
