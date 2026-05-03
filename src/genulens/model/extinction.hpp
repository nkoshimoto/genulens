#pragma once

namespace genulens::model {

struct BandExtinction {
    double i_band = 0.0;
    double k_band = 0.0;
    double color_vi = 0.0;
};

class ExponentialDustExtinction {
public:
    ExponentialDustExtinction(double scale_pc, double reference_distance_pc,
                              double ai_reference, double ak_reference, double evi_reference);

    double scale_pc() const { return scale_pc_; }
    double reference_distance_pc() const { return reference_distance_pc_; }
    double ai0() const { return ai0_; }
    double ak0() const { return ak0_; }
    double evi0() const { return evi0_; }
    BandExtinction at_distance(double distance_pc) const;
    double distance_modulus_term(double distance_pc) const;

private:
    double scale_pc_ = 0.0;
    double reference_distance_pc_ = 0.0;
    double ai0_ = 0.0;
    double ak0_ = 0.0;
    double evi0_ = 0.0;
};

double source_distance_weight(double distance_pc, double gamma_ds);

} // namespace genulens::model
