#include "genulens/simulation/scientific_cli.hpp"

#include "genulens/simulation/scientific_engine.hpp"

namespace genulens {

int run_scientific_cli(int argc, char **argv)
{
    ScientificEngine engine;
    return engine.run(argc, argv);
}

} // namespace genulens
