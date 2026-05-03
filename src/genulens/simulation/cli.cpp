#include "genulens/simulation/cli.hpp"

#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/sampler.hpp"

namespace genulens {

int run_cli(int argc, char **argv)
{
    Initializer initializer;
    auto context = initializer.create_context();
    Sampler sampler;
    return sampler.run_cli(context, argc, argv);
}

} // namespace genulens
