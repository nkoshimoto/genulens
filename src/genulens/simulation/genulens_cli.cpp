#include "genulens/simulation/genulens_cli.hpp"

#include "genulens/simulation/genulens_initializer.hpp"
#include "genulens/simulation/genulens_sampler.hpp"

namespace genulens {

int run_genulens_cli(int argc, char **argv)
{
    GenulensInitializer initializer;
    auto context = initializer.create_context();
    GenulensSampler sampler;
    return sampler.run_cli(context, argc, argv);
}

} // namespace genulens
