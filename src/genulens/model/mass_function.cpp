#include "genulens/model/mass_function.hpp"

#include "genulens/model/parameters.hpp"

#include <cmath>

namespace genulens::model {

double broken_power_law_imf(double mass_msun, const IMFParameters &params)
{
    if (mass_msun <= 0.0) return 0.0;
    if (mass_msun < params.m3) return std::pow(mass_msun, params.alpha4);
    if (mass_msun < params.m2) return std::pow(mass_msun, params.alpha3);
    if (mass_msun < params.m1) return std::pow(mass_msun, params.alpha2);
    if (mass_msun < params.m0) return std::pow(mass_msun, params.alpha1);
    return std::pow(mass_msun, params.alpha0);
}

double broken_power_law_imf(double mass_msun)
{
    return broken_power_law_imf(mass_msun, default_model_parameters().imf);
}

} // namespace genulens::model
