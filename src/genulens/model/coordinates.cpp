#include "genulens/model/coordinates.hpp"

#include "genulens/constants.hpp"

#include <cmath>

namespace genulens::model {

namespace {

void cross(double *c, const double *a, const double *b)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double dot(const double *a, const double *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void norm_vec(double *a)
{
    const double norm = std::sqrt(dot(a, a));
    if (norm == 0.0) return;
    a[0] /= norm;
    a[1] /= norm;
    a[2] /= norm;
}

} // namespace

CoordinateTransformer::CoordinateTransformer(std::array<double, 3> sgr_a_offset)
    : sgr_a_offset_(sgr_a_offset)
{
}

std::array<double, 3> CoordinateTransformer::distance_l_b_to_xyz(double distance_pc, double l_deg, double b_deg,
                                                                 double sun_radius_pc) const
{
    const double cosbsun = std::cos(kSunZPc / sun_radius_pc);
    const double sinbsun = std::sin(kSunZPc / sun_radius_pc);
    const double b = b_deg * kPi / 180.0;
    const double l = l_deg * kPi / 180.0;
    const double xtmp = sun_radius_pc - distance_pc * std::cos(b) * std::cos(l);
    const double ytmp = distance_pc * std::cos(b) * std::sin(l);
    const double ztmp = distance_pc * std::sin(b);
    return {
        xtmp - sgr_a_offset_[0],
        ytmp - sgr_a_offset_[1],
        ztmp * cosbsun + xtmp * sinbsun - sgr_a_offset_[2],
    };
}

PositionAngle CoordinateTransformer::position_angle(double l_deg, double b_deg)
{
    const double l_np = 122.9320;
    const double b_np = 27.1284;
    const double l = l_deg * kPi / 180.0;
    const double b = b_deg * kPi / 180.0;
    const double lnp = l_np * kPi / 180.0;
    const double bnp = b_np * kPi / 180.0;

    double nvector[3] = {std::cos(b) * std::cos(l), std::cos(b) * std::sin(l), std::sin(b)};
    double gal_np[3] = {0.0, 0.0, 1.0};
    double eq_np[3] = {std::cos(bnp) * std::cos(lnp), std::cos(bnp) * std::sin(lnp), std::sin(bnp)};

    double elvector[3] = {};
    double ebvector[3] = {};
    double eevector[3] = {};
    double envector[3] = {};
    cross(elvector, gal_np, nvector);
    norm_vec(elvector);
    cross(ebvector, nvector, elvector);
    cross(eevector, eq_np, nvector);
    norm_vec(eevector);
    cross(envector, nvector, eevector);

    double crosstmp[3] = {};
    PositionAngle pa;
    pa.cos_pa = dot(elvector, eevector);
    cross(crosstmp, elvector, eevector);
    pa.sin_pa = -dot(nvector, crosstmp);
    pa.degrees = 180.0 * std::atan2(pa.sin_pa, pa.cos_pa) / kPi;
    return pa;
}

} // namespace genulens::model

