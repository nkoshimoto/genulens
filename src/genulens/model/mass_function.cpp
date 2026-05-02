#include "genulens/model/mass_function.hpp"

#include <cmath>

namespace genulens::model {

double broken_power_law_imf(double mass_msun)
{
    if (mass_msun <= 0.0) return 0.0;
    if (mass_msun < 0.01) return std::pow(mass_msun, -0.18);
    if (mass_msun < 0.08) return std::pow(mass_msun, -0.18);
    if (mass_msun < 0.86) return std::pow(mass_msun, -1.13);
    return std::pow(mass_msun, -2.32);
}

} // namespace genulens::model

