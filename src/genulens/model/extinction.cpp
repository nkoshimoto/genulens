#include "genulens/model/extinction.hpp"

#include <cmath>

namespace genulens::model {

ExponentialDustExtinction::ExponentialDustExtinction(double scale_pc, double reference_distance_pc,
                                                     const ReferenceBandExtinction &reference)
    : scale_pc_(scale_pc),
      reference_distance_pc_(reference_distance_pc)
{
    const double denominator = (scale_pc_ > 0.0 && reference_distance_pc_ > 0.0)
                                   ? 1.0 - std::exp(-reference_distance_pc_ / scale_pc_)
                                   : 0.0;
    av0_ = denominator != 0.0 ? reference.v_band / denominator : 0.0;
    ai0_ = denominator != 0.0 ? reference.i_band / denominator : 0.0;
    aj0_ = denominator != 0.0 ? reference.j_band / denominator : 0.0;
    ah0_ = denominator != 0.0 ? reference.h_band / denominator : 0.0;
    ak0_ = denominator != 0.0 ? reference.k_band / denominator : 0.0;
    evi0_ = denominator != 0.0 ? reference.color_vi / denominator : 0.0;
    f0870_ = denominator != 0.0 ? reference.f087_band / denominator : 0.0;
    f1460_ = denominator != 0.0 ? reference.f146_band / denominator : 0.0;
    f2130_ = denominator != 0.0 ? reference.f213_band / denominator : 0.0;
}

BandExtinction ExponentialDustExtinction::at_distance(double distance_pc) const
{
    const double factor = 1.0 - std::exp(-distance_pc / scale_pc_);
    return {av0_ * factor, ai0_ * factor, aj0_ * factor, ah0_ * factor,
            ak0_ * factor, evi0_ * factor, f0870_ * factor, f1460_ * factor,
            f2130_ * factor};
}

double ExponentialDustExtinction::distance_modulus_term(double distance_pc) const
{
    return 5.0 * std::log10(0.1 * distance_pc);
}

double source_distance_weight(double distance_pc, double gamma_ds)
{
    return distance_pc > 0.0 ? std::pow(distance_pc / 8000.0, gamma_ds) : 0.0;
}

namespace {

double alam_av_wang_chen_2019(double wavelength_nm)
{
    if (wavelength_nm < 1000.0) {
        const double y1 = 1000.0 / wavelength_nm - 1.82;
        const double coeffs[7] = {0.7499, -0.1086, -0.08909, 0.02905,
                                  0.01069, 0.001707, -0.001002};
        double yi = 1.0;
        double alam_av = 1.0;
        for (double coeff : coeffs) {
            yi *= y1;
            alam_av += coeff * yi;
        }
        return alam_av;
    }
    return 0.3722 * std::pow(1000.0 / wavelength_nm, 2.07);
}

double genstars_alam_per_ejk(double l_deg, double b_deg, double wavelength_nm,
                             int extinction_law)
{
    if (extinction_law == 2) {
        const double ejk_av = alam_av_wang_chen_2019(1254.0) -
                              alam_av_wang_chen_2019(2149.0);
        return alam_av_wang_chen_2019(wavelength_nm) / ejk_av;
    }

    const double ejk_ai = (l_deg > 0.0 && b_deg > 0.0) ? 3.65
                         : (l_deg < 0.0 && b_deg > 0.0) ? 3.77
                         : (l_deg > 0.0 && b_deg < 0.0) ? 3.97
                         : (l_deg < 0.0 && b_deg < 0.0) ? 3.82
                         : 3.86;
    const double ai_av = (l_deg > 0.0 && b_deg > 0.0) ? 1.80
                         : (l_deg < 0.0 && b_deg > 0.0) ? 1.81
                         : (l_deg > 0.0 && b_deg < 0.0) ? 1.82
                         : (l_deg < 0.0 && b_deg < 0.0) ? 1.82
                         : 1.82;
    const double alpha_i_to_v = (l_deg > 0.0 && b_deg > 0.0) ? 1.54
                                : (l_deg < 0.0 && b_deg > 0.0) ? 1.55
                                : (l_deg > 0.0 && b_deg < 0.0) ? 1.57
                                : (l_deg < 0.0 && b_deg < 0.0) ? 1.56
                                : 1.56;
    const double ejk_av = ejk_ai * ai_av;

    double ejk_ak = 0.0;
    double ehk_ak = 0.0;
    double alpha_j_to_i = 0.0;
    double alpha_ir = 0.0;
    double f_to_ejk_vvv = 1.0;
    double wavelength0[5] = {549.056, 805.988, 0.0, 0.0, 0.0};
    if (extinction_law == 1) {
        ejk_ak = (l_deg > 0.0 && b_deg > 0.0) ? 0.497
                 : (l_deg < 0.0 && b_deg > 0.0) ? 0.494
                 : (l_deg > 0.0 && b_deg < 0.0) ? 0.534
                 : (l_deg < 0.0 && b_deg < 0.0) ? 0.587
                 : 0.528;
        ehk_ak = (l_deg > 0.0 && b_deg > 0.0) ? 1.64
                 : (l_deg < 0.0 && b_deg > 0.0) ? 1.48
                 : (l_deg > 0.0 && b_deg < 0.0) ? 1.54
                 : (l_deg < 0.0 && b_deg < 0.0) ? 1.63
                 : 1.61;
        alpha_j_to_i = (l_deg > 0.0 && b_deg > 0.0) ? 2.07
                       : (l_deg < 0.0 && b_deg > 0.0) ? 2.15
                       : (l_deg > 0.0 && b_deg < 0.0) ? 2.21
                       : (l_deg < 0.0 && b_deg < 0.0) ? 2.05
                       : 2.12;
        alpha_ir = 2.0;
        wavelength0[2] = 1240.0;
        wavelength0[3] = 1664.0;
        wavelength0[4] = 2164.0;
        f_to_ejk_vvv = 0.970;
    } else {
        ejk_ak = (l_deg > 0.0 && b_deg > 0.0) ? 0.390
                 : (l_deg < 0.0 && b_deg > 0.0) ? 0.384
                 : (l_deg > 0.0 && b_deg < 0.0) ? 0.464
                 : (l_deg < 0.0 && b_deg < 0.0) ? 0.415
                 : 0.428;
        ehk_ak = (l_deg > 0.0 && b_deg > 0.0) ? 1.02
                 : (l_deg < 0.0 && b_deg > 0.0) ? 0.97
                 : (l_deg > 0.0 && b_deg < 0.0) ? 1.30
                 : (l_deg < 0.0 && b_deg < 0.0) ? 1.21
                 : 1.10;
        alpha_j_to_i = (l_deg > 0.0 && b_deg > 0.0) ? 2.18
                       : (l_deg < 0.0 && b_deg > 0.0) ? 2.26
                       : (l_deg > 0.0 && b_deg < 0.0) ? 2.26
                       : (l_deg < 0.0 && b_deg < 0.0) ? 2.25
                       : 2.25;
        alpha_ir = 2.47;
        wavelength0[2] = 1254.0;
        wavelength0[3] = 1646.0;
        wavelength0[4] = 2149.0;
        f_to_ejk_vvv = 1.0;
    }

    const double ejk_aj = ejk_ak + 1.0;
    const double ejk_ah = ejk_ak * (1.0 / ehk_ak + 1.0);
    const double alam0[5] = {ejk_av, ejk_ai, ejk_aj / f_to_ejk_vvv,
                             ejk_ah / f_to_ejk_vvv, ejk_ak / f_to_ejk_vvv};
    int closest = 0;
    double closest_distance = std::abs(wavelength_nm - wavelength0[0]);
    for (int i = 1; i < 5; ++i) {
        const double distance = std::abs(wavelength_nm - wavelength0[i]);
        if (distance < closest_distance) {
            closest = i;
            closest_distance = distance;
        }
    }
    const double alpha = (wavelength_nm < wavelength0[1]) ? alpha_i_to_v
                         : (wavelength_nm < 1000.0) ? alpha_j_to_i
                         : alpha_ir;
    return std::pow(wavelength0[closest] / wavelength_nm, alpha) * alam0[closest];
}

} // namespace

ReferenceBandExtinction genstars_reference_extinction(double l_deg, double b_deg,
                                                      double ejk_reference,
                                                      int extinction_law)
{
    ReferenceBandExtinction ref;
    ref.v_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 549.056, extinction_law);
    ref.i_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 805.988, extinction_law);
    ref.j_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 1240.0, extinction_law);
    ref.h_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 1664.0, extinction_law);
    ref.k_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 2164.0, extinction_law);
    ref.color_vi = ref.v_band - ref.i_band;
    ref.f087_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 867.590, extinction_law);
    ref.f146_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 1367.793, extinction_law);
    ref.f213_band = ejk_reference * genstars_alam_per_ejk(l_deg, b_deg, 2112.465, extinction_law);
    return ref;
}

} // namespace genulens::model
