#pragma once

#include <vector>
#include "genulens/model/extinction.hpp"
#include "genulens/simulation/run_context.hpp"

namespace genulens {

struct LineOfSightDensityGridConfig {
    int Dmax = 16000;
    double l = 0.0, b = 0.0;

    int BHhd = 0, BHhb = 0, fixRhdBH = 0;
    double vkickBH = 100.0, RhdBH0 = 9660.0, betaBH = 0.820;
    int UseSigBH = 0;

    double AI0 = 0.0, EVI0 = 0.0, gammaDs = 0.5;
    double Isst = 14.0, Isen = 21.0, VIsst = 0.0, VIsen = 0.0;
    double wtD_L = 0.0;
    bool check_D = false;

    int npri = 10;
    bool printBHfac = false;
    bool printrhoS = false;

    const model::ExponentialDustExtinction *extinction = nullptr;
};

class LineOfSightDensityGrid {
public:
    void build(RunContext &ctx, const LineOfSightDensityGridConfig &cfg);

    int sample_source_component(double ran) const;
    int sample_lens_component(int nbinDlmin, int nbinDlmax, double ran) const;
    double sample_source_distance(int i_s, double ran) const;
    double sample_lens_distance(int i_l, int nbinDlmin, int nbinDlmax, double ran) const;

    void get_lens_importance_bins(double D_s,
                                   double thetaEmin, double thetaEmax,
                                   double piEmin, double piEmax,
                                   int &nbinDlmin, int &nbinDlmax) const;
    double importance_weight(int nbinDlmin, int nbinDlmax) const;
    double fBH_at_distance(int component, double D) const;

    int nbin() const { return nbin_; }
    int ncomp() const { return ncomp_; }
    double dD() const { return dD_; }
    double total_source_count() const { return nallS_; }

    const double *D_data() const { return D_.data(); }
    const double *cumu_rho_all_S_data() const { return cumu_rho_all_S_.data(); }
    const double *cumu_rho_all_L_data() const { return cumu_rho_all_L_.data(); }
    const double *cumu_rho_S_data(int i) const { return cumu_rho_S_[i].data(); }
    const double *cumu_rho_L_data(int i) const { return cumu_rho_L_[i].data(); }
    const double *rhoD_S_data(int i) const { return rhoD_S_[i].data(); }
    const double *rhoD_L_data(int i) const { return rhoD_L_[i].data(); }

private:
    int nbin_ = 0, ncomp_ = 0;
    double dD_ = 0.0;
    double nallS_ = 0.0;

    std::vector<double> D_;
    std::vector<double> cumu_rho_all_S_, cumu_rho_all_L_;
    std::vector<std::vector<double>> rhoD_S_, rhoD_L_;
    std::vector<std::vector<double>> cumu_rho_S_, cumu_rho_L_;
    std::vector<std::vector<double>> fBH_;
    std::vector<std::vector<int>> ibinptiles_S_, ibinptiles_L_;
};

} // namespace genulens
