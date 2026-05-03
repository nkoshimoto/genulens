#include "genulens/simulation/legacy_cli.hpp"

int genulens_legacy_main(int argc, char **argv);

namespace genulens {

int run_legacy_cli(int argc, char **argv)
{
    return genulens_legacy_main(argc, argv);
}

} // namespace genulens

