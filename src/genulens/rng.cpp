#include "genulens/rng.hpp"

#include <gsl/gsl_randist.h>

#include <stdexcept>
#include <utility>

namespace genulens {

RandomEngine::RandomEngine(unsigned long seed)
{
    gsl_rng_env_setup();
    rng_ = gsl_rng_alloc(gsl_rng_default);
    if (!rng_) {
        throw std::runtime_error("failed to allocate GSL RNG");
    }
    gsl_rng_set(rng_, seed);
}

RandomEngine::~RandomEngine()
{
    gsl_rng_free(rng_);
}

RandomEngine::RandomEngine(RandomEngine &&other) noexcept
    : rng_(std::exchange(other.rng_, nullptr))
{
}

RandomEngine &RandomEngine::operator=(RandomEngine &&other) noexcept
{
    if (this != &other) {
        gsl_rng_free(rng_);
        rng_ = std::exchange(other.rng_, nullptr);
    }
    return *this;
}

double RandomEngine::uniform()
{
    return gsl_rng_uniform(rng_);
}

double RandomEngine::gaussian()
{
    return gsl_ran_ugaussian(rng_);
}

} // namespace genulens

