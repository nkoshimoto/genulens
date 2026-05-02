#include "genulens/model/galactic_model.hpp"
#include "genulens/options.hpp"
#include "genulens/tools/rho_profile.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    const auto cfg = genulens::config_from_cli(argc, argv);
    const genulens::GalacticModel model(cfg);
    for (const auto &row : genulens::tools::rho_profile(model, 100.0, 16000.0, 100.0)) {
        std::cout << row.distance_pc << '\t' << row.density.total() << '\n';
    }
    return 0;
}

