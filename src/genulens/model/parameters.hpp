#pragma once

namespace genulens::model {

struct IMFParameters {
    double m0 = 1.0;
    double m1 = 0.859770466578045;
    double m2 = 0.08;
    double m3 = 0.01;
    double ml = 0.001;
    double mu = 120.0;
    double alpha0 = -2.32279457078378;
    double alpha1 = -2.32279457078378;
    double alpha2 = -1.13449983242887;
    double alpha3 = -0.175862190587576;
    double alpha4 = -0.175862190587576;
};

struct ModelParameters {
    IMFParameters imf;
};

const ModelParameters &default_model_parameters();

} // namespace genulens::model

