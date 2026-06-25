#pragma once

namespace genulens::model {

struct BandExtinction {
    double v_band = 0.0;
    double i_band = 0.0;
    double j_band = 0.0;
    double h_band = 0.0;
    double k_band = 0.0;
    double color_vi = 0.0;
    double f087_band = 0.0;
    double f146_band = 0.0;
    double f213_band = 0.0;
};

struct ReferenceBandExtinction {
    double v_band = 0.0;
    double i_band = 0.0;
    double j_band = 0.0;
    double h_band = 0.0;
    double k_band = 0.0;
    double color_vi = 0.0;
    double f087_band = 0.0;
    double f146_band = 0.0;
    double f213_band = 0.0;
};

class ExponentialDustExtinction {
public:
    ExponentialDustExtinction(double scale_pc, double reference_distance_pc,
                              const ReferenceBandExtinction &reference);

    double scale_pc() const { return scale_pc_; }
    double reference_distance_pc() const { return reference_distance_pc_; }
    double av0() const { return av0_; }
    double ai0() const { return ai0_; }
    double aj0() const { return aj0_; }
    double ah0() const { return ah0_; }
    double ak0() const { return ak0_; }
    double evi0() const { return evi0_; }
    double f0870() const { return f0870_; }
    double f1460() const { return f1460_; }
    double f2130() const { return f2130_; }
    BandExtinction at_distance(double distance_pc) const;
    double distance_modulus_term(double distance_pc) const;

private:
    double scale_pc_ = 0.0;
    double reference_distance_pc_ = 0.0;
    double av0_ = 0.0;
    double ai0_ = 0.0;
    double aj0_ = 0.0;
    double ah0_ = 0.0;
    double ak0_ = 0.0;
    double evi0_ = 0.0;
    double f0870_ = 0.0;
    double f1460_ = 0.0;
    double f2130_ = 0.0;
};

ReferenceBandExtinction genstars_reference_extinction(double l_deg, double b_deg,
                                                      double ejk_reference,
                                                      int extinction_law);
double source_distance_weight(double distance_pc, double gamma_ds);

} // namespace genulens::model
