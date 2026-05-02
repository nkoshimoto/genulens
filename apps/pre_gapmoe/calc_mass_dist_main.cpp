#include "genulens/model/galactic_model.hpp"
#include "genulens/options.hpp"
#include "genulens/tools/mass_dist.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    const auto cfg = genulens::config_from_cli(argc, argv);
    const genulens::GalacticModel model(cfg);
    for (const auto &row : genulens::tools::mass_distribution(model, 1000)) {
        std::cout << row.log_mass << '\t' << row.dndlogm << '\n';
    }
    return 0;
}

