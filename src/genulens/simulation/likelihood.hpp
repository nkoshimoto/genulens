#pragma once

#include "genulens/types.hpp"

#include <functional>

namespace genulens {

using LikelihoodFunction = std::function<double(const Event &)>;

class GaussianLikelihood {
public:
    GaussianLikelihood(double observed_tE, double sigma_tE);
    double operator()(const Event &event) const;

private:
    double observed_tE_;
    double sigma_tE_;
};

LikelihoodFunction default_likelihood(double observed_tE, double sigma_tE);

} // namespace genulens

