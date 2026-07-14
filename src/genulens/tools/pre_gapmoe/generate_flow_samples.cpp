#include "genulens/cli/option.h"
#include "genulens/simulation/initialize.hpp"
#include "genulens/simulation/sampler.hpp"

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool has_option(int argc, char **argv, const char *name)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], name) == 0) return true;
    }
    return false;
}

bool wants_help(int argc, char **argv)
{
    return has_option(argc, argv, "--help") || has_option(argc, argv, "-h") ||
           has_option(argc, argv, "help");
}

void print_usage(const char *program)
{
    std::cout
        << "Usage: " << program << " SOURCEGROUP <3|4> [genulens options]\n\n"
        << "Generate full weighted microlensing events for conditional-flow training.\n"
        << "The source is restricted to NSD (3) or stellar halo (4), while lens\n"
        << "components retain the standard Galactic mixture. Output uses the normal\n"
        << "genulens VERBOSITY=3 event schema.\n\n"
        << "Examples:\n"
        << "  " << program << " SOURCEGROUP 4 l 1 b -3.9 Nsimu 100000 seed 42\n"
        << "  " << program << " SOURCEGROUP 3 NSD 1 l 0 b 0 Nsimu 100000 seed 43\n";
}

} // namespace

int main(int argc, char **argv)
{
    if (wants_help(argc, argv)) {
        print_usage(argv[0]);
        return 0;
    }

    const int source_group = getOptioni(argc, argv, "SOURCEGROUP", 1, -1);
    if (source_group != 3 && source_group != 4) {
        std::cerr << "SOURCEGROUP must be 3 (NSD) or 4 (stellar halo)\n";
        return 2;
    }
    if (has_option(argc, argv, "VERBOSITY") &&
        getOptioni(argc, argv, "VERBOSITY", 1, 3) != 3) {
        std::cerr << "generate_flow_samples requires VERBOSITY 3\n";
        return 2;
    }

    const int source_component = source_group == 3 ? 9 : 10;
    std::vector<std::string> args;
    args.emplace_back(argv[0]);
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "SOURCEGROUP") == 0) {
            ++i;
            continue;
        }
        args.emplace_back(argv[i]);
    }
    if (!has_option(argc, argv, "VERBOSITY")) {
        args.emplace_back("VERBOSITY");
        args.emplace_back("3");
    }
    if (source_group == 3 && !has_option(argc, argv, "NSD")) {
        args.emplace_back("NSD");
        args.emplace_back("1");
    }
    args.emplace_back("__PRE_GAPMOE_SOURCE_COMPONENT");
    args.emplace_back(std::to_string(source_component));

    std::vector<char *> forwarded;
    forwarded.reserve(args.size());
    for (auto &arg : args) forwarded.push_back(arg.data());

    try {
        // Keep the Python pre_gapmoe table parser self-describing.  The
        // standard CLI header is human-readable but does not use its
        // ``# Columns:`` convention.
        std::cout
            << "# Columns:\n"
            << "# weight ML DL DS tE thetaE piE piEN piEE mu_rel muSl muSb "
               "IL KL source_component lens_component remnant\n";
        genulens::Initializer initializer;
        auto context = initializer.create_context();
        initializer.initialize_rng(
            context, static_cast<int>(forwarded.size()), forwarded.data());
        return genulens::Sampler().run_cli(
            context, static_cast<int>(forwarded.size()), forwarded.data());
    } catch (const std::exception &error) {
        std::cerr << "generate_flow_samples: " << error.what() << '\n';
        return 3;
    }
}
