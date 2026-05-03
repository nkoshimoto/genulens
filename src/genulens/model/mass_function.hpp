#pragma once

namespace genulens::model {

struct IMFParameters;

double broken_power_law_imf(double mass_msun, const IMFParameters &params);
double broken_power_law_imf(double mass_msun);

} // namespace genulens::model
