#pragma once

namespace genulens::math {

class Interpolation {
public:
    static double linear(int n, double *x, double *y, double xin);
    static double linear_from_index(int n, double *x, double *y, double xin, int *ist);
    static double linear_with_upper_index(int n, double *x, double *y, double xin, int *khi);
    static double uniform_grid(int n, double *values, double x_start, double dx, double x_req);
    static double bilinear(int nx, int ny, double **values, double x_start, double y_start,
                           double dx, double dy, double x_req, double y_req);
    static void bilinear_coefficients(int nx, int ny, double *coefficients,
                                      double x_start, double y_start, double dx, double dy,
                                      double x_req, double y_req);
    static double inverse_cumulative_linear_density(int n, double *x, double *cumulative,
                                                    double *density, double requested,
                                                    int start_index, int inverse_direction);
};

} // namespace genulens::math

