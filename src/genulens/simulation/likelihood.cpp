#include "genulens/simulation/likelihood.hpp"

#include <cmath>

namespace genulens {

GaussianLikelihood::GaussianLikelihood(double observed_tE, double sigma_tE)
    : observed_tE_(observed_tE), sigma_tE_(sigma_tE)
{
}

double GaussianLikelihood::operator()(const Event &event) const
{
    if (sigma_tE_ <= 0.0) return 1.0;
    const double pull = (event.tE - observed_tE_) / sigma_tE_;
    return std::exp(-0.5 * pull * pull);
}

LikelihoodFunction default_likelihood(double observed_tE, double sigma_tE)
{
    return GaussianLikelihood(observed_tE, sigma_tE);
}

} // namespace genulens

