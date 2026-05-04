#include "genulens/simulation/velocity_distribution.hpp"

#include <array>

#include "genulens/simulation/internal/runtime.hpp"

namespace genulens {

std::array<double, 3> VelocityDistribution::sample(
        int component, double age, double D, double l, double b) {
    void get_vxyz_ran(RunContext &ctx, double *vxyz, int i, double tau, double D, double lD, double bD);
    double vxyz[3] = {};
    get_vxyz_ran(*ctx_, vxyz, component, age, D, l, b);
    return {vxyz[0], vxyz[1], vxyz[2]};
}

} // namespace genulens
