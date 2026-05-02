#pragma once

#include <gsl/gsl_rng.h>

namespace genulens {

class RandomEngine {
public:
    explicit RandomEngine(unsigned long seed = 12304357UL);
    ~RandomEngine();

    RandomEngine(const RandomEngine &) = delete;
    RandomEngine &operator=(const RandomEngine &) = delete;

    RandomEngine(RandomEngine &&other) noexcept;
    RandomEngine &operator=(RandomEngine &&other) noexcept;

    double uniform();
    double gaussian();
    gsl_rng *native() { return rng_; }
    const gsl_rng *native() const { return rng_; }

private:
    gsl_rng *rng_ = nullptr;
};

} // namespace genulens

