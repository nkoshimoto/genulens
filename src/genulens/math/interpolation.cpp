#include "genulens/math/interpolation.hpp"

#include <cmath>

namespace genulens::math {

double Interpolation::linear(int n, double *x, double *y, double xin)
{
    double xmin, xmax;
    if (x[0] < x[n - 1]) {
        xmin = x[0];
        xmax = x[n - 1];
    } else {
        xmin = x[n - 1];
        xmax = x[0];
    }
    if (xmin > xin || xmax < xin) return 0.0;

    for (int i = 1; i < n; ++i) {
        if ((x[i] <= xin && x[i - 1] >= xin) || (x[i] >= xin && x[i - 1] <= xin)) {
            return (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (xin - x[i - 1]) + y[i - 1];
        }
    }
    return 0.0;
}

double Interpolation::linear_from_index(int n, double *x, double *y, double xin, int *ist)
{
    for (int i = *ist; i < n; ++i) {
        if (i == 0) continue;
        if ((x[i] <= xin && x[i - 1] >= xin) || (x[i] >= xin && x[i - 1] <= xin)) {
            *ist = i;
            return (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (xin - x[i - 1]) + y[i - 1];
        }
    }
    return 0.0;
}

double Interpolation::linear_with_upper_index(int n, double *x, double *y, double xin, int *khi)
{
    double xmin, xmax;
    if (x[0] < x[n - 1]) {
        xmin = x[0];
        xmax = x[n - 1];
    } else {
        xmin = x[n - 1];
        xmax = x[0];
    }
    if (xmin > xin || xmax < xin) return 0.0;

    int klo;
    if (*khi > 0) {
        klo = *khi - 1;
    } else {
        klo = 0;
        *khi = n - 1;
        while (*khi - klo > 1) {
            int k = (*khi + klo) >> 1;
            if (x[k] > xin) *khi = k;
            else klo = k;
        }
    }
    const double h = x[*khi] - x[klo];
    if (h == 0.0) return 0.0;
    const double a = (x[*khi] - xin) / h;
    const double b = (xin - x[klo]) / h;
    return a * y[klo] + b * y[*khi];
}

double Interpolation::uniform_grid(int n, double *values, double x_start, double dx, double x_req)
{
    const int ix = static_cast<int>((x_req - x_start) / dx);
    const double xres = (x_req - x_start) / dx - ix;
    if (ix < 0 || ix > n - 1) return 0.0;
    if (ix + 1 > n - 1) return values[ix];
    return values[ix] * (1.0 - xres) + values[ix + 1] * xres;
}

double Interpolation::bilinear(int nx, int ny, double **values, double x_start, double y_start,
                               double dx, double dy, double x_req, double y_req)
{
    const int ix = static_cast<int>((x_req - x_start) / dx);
    const double xres = (x_req - x_start) / dx - ix;
    const int iy = static_cast<int>((y_req - y_start) / dy);
    const double yres = (y_req - y_start) / dy - iy;
    if (ix < 0 || ix > nx - 1 || iy < 0 || iy > ny - 1) return 0.0;
    if (ix + 1 > nx - 1 && iy + 1 > ny - 1) return values[ix][iy];
    if (ix + 1 > nx - 1) return values[ix][iy] * (1.0 - yres) + values[ix][iy + 1] * yres;
    if (iy + 1 > ny - 1) return values[ix][iy] * (1.0 - xres) + values[ix + 1][iy] * xres;
    return (1.0 - xres) * (1.0 - yres) * values[ix][iy] +
           xres * (1.0 - yres) * values[ix + 1][iy] +
           (1.0 - xres) * yres * values[ix][iy + 1] +
           xres * yres * values[ix + 1][iy + 1];
}

void Interpolation::bilinear_coefficients(int nx, int ny, double *coefficients,
                                          double x_start, double y_start, double dx, double dy,
                                          double x_req, double y_req)
{
    const int ix = static_cast<int>((x_req - x_start) / dx);
    const double xres = (x_req - x_start) / dx - ix;
    const int iy = static_cast<int>((y_req - y_start) / dy);
    const double yres = (y_req - y_start) / dy - iy;
    if (ix < 0 || ix > nx - 1 || iy < 0 || iy > ny - 1) {
        coefficients[0] = coefficients[1] = coefficients[2] = coefficients[3] = 0.0;
    } else if (ix + 1 > nx - 1 && iy + 1 > ny - 1) {
        coefficients[0] = 1.0;
        coefficients[1] = coefficients[2] = coefficients[3] = 0.0;
    } else if (ix + 1 > nx - 1) {
        coefficients[0] = 1.0 - yres;
        coefficients[2] = yres;
        coefficients[1] = coefficients[3] = 0.0;
    } else if (iy + 1 > ny - 1) {
        coefficients[0] = 1.0 - xres;
        coefficients[1] = xres;
        coefficients[2] = coefficients[3] = 0.0;
    } else {
        coefficients[0] = (1.0 - xres) * (1.0 - yres);
        coefficients[1] = xres * (1.0 - yres);
        coefficients[2] = (1.0 - xres) * yres;
        coefficients[3] = xres * yres;
    }
}

double Interpolation::inverse_cumulative_linear_density(int n, double *x, double *cumulative,
                                                        double *density, double requested,
                                                        int start_index, int inverse_direction)
{
    if (cumulative[0] > requested || cumulative[n - 1] < requested) return 0.0;
    if (start_index < 1) start_index = 1;

    const auto solve = [&](int i) {
        const double a = 0.5 * (density[i] - density[i - 1]) / (x[i] - x[i - 1]);
        const double b = density[i - 1] - 2.0 * a * x[i - 1];
        const double c = a * x[i - 1] * x[i - 1] - density[i - 1] * x[i - 1] + cumulative[i - 1] - requested;
        return (a != 0.0) ? (-b + std::sqrt(b * b - 4.0 * a * c)) * 0.5 / a
                          : (x[i] - x[i - 1]) / (cumulative[i] - cumulative[i - 1]) * (requested - cumulative[i - 1]) + x[i - 1];
    };

    if (inverse_direction == 0) {
        for (int i = start_index; i < n; ++i) {
            if ((cumulative[i] <= requested && cumulative[i - 1] > requested) ||
                (cumulative[i] >= requested && cumulative[i - 1] < requested)) {
                return solve(i);
            }
        }
    } else {
        for (int i = start_index; i > 0; --i) {
            if ((cumulative[i] <= requested && cumulative[i - 1] > requested) ||
                (cumulative[i] >= requested && cumulative[i - 1] < requested)) {
                return solve(i);
            }
        }
    }
    return 0.0;
}

} // namespace genulens::math

