#include "genulens/io/output.hpp"
#include "genulens/options.hpp"
#include "genulens/simulation/simulator.hpp"

#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help" || arg == "help") {
            std::cout << "Usage: genulens_core_cli [l <deg>] [b <deg>] [Nsimu <n>] [seed <n>]\n";
            return 0;
        }
    }
    const auto cfg = genulens::config_from_cli(argc, argv);
    genulens::write_events_tsv(std::cout, genulens::simulate(cfg));
    return 0;
}

