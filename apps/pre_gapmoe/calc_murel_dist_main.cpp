#include "genulens/rng.hpp"
#include "genulens/tools/murel_dist.hpp"

#include <iostream>

int main()
{
    genulens::RandomEngine rng;
    for (const auto &row : genulens::tools::murel_distribution(rng, 10000, 300.0, 0.5)) {
        std::cout << row.murel << '\t' << row.probability << '\n';
    }
    return 0;
}

