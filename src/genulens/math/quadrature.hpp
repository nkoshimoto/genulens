#pragma once

namespace genulens::math {

class NewtonCotes {
public:
    static int coefficients(int requested_order, double *locations, double *weights);
};

} // namespace genulens::math
