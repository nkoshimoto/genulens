#include "genulens/model/coordinates.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"

namespace gmodel = genulens::model;
#include "genulens/simulation/internal/runtime.hpp"
namespace genulens {

void cross(double *c, double *a, double *b){
   c[0] = a[1]*b[2] - b[1]*a[2];
   c[1] = a[2]*b[0] - b[2]*a[0];
   c[2] = a[0]*b[1] - b[0]*a[1];
}
//-----get c= a . b------------------
double dot(double *a, double *b){
   return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
//---- normalize a vector ----
void norm_vec(double *a){
  double norm = dot(a, a);
  for (int i=0;i<3;i++){
    a[i] = a[i]/sqrt(norm);
  }
}
void Dlb2xyz(double D, double lD, double bD, double Rsun, double *xyz)
/*------------------------------------------------------------*/
/*  Give (x,y,z) for a given D, lD, bD  */
/*------------------------------------------------------------*/
{
  const std::array<double, 3> offset = {xyzSgrA[0], xyzSgrA[1], xyzSgrA[2]};
  const auto value = gmodel::CoordinateTransformer(offset).distance_l_b_to_xyz(D, lD, bD, Rsun);
  xyz[0] = value[0];
  xyz[1] = value[1];
  xyz[2] = value[2];
}


/*----------------------------------------------------------------*/
/*                      for general use                           */
/*----------------------------------------------------------------*/
//---- get parameters for trapezoidal integral -------
// Values from http://midarekazu.g2.xrea.com/Newton-Cotes.html
int get_p_integral(int nji, double *ls, double *ks)
{
  return genulens::math::NewtonCotes::coefficients(nji, ls, ks);
}
/*----------------------------------------------------------------*/
//---- getx2y for linear interpolation
double getx2y(int n, double *x, double *y, double xin)
{
   return genulens::math::Interpolation::linear(n, x, y, xin);
}
//---- getx2y_ist for linear interpolation
double getx2y_ist(int n, double *x, double *y, double xin, int *ist)
{
   return genulens::math::Interpolation::linear_from_index(n, x, y, xin, ist);
}
//---- getx2y_khi for linear interpolation
double getx2y_khi(int n, double *x, double *y, double xin, int *khi)
{
   return genulens::math::Interpolation::linear_with_upper_index(n, x, y, xin, khi);
}
//---------------
double interp_x(int n, double *F, double xst, double dx, double xreq) // just for this code
{
  return genulens::math::Interpolation::uniform_grid(n, F, xst, dx, xreq);
}
//---------------
double interp_xy(int nx, int ny, double **F, double xst, double yst, double dx, double dy, double xreq, double yreq) // just for this code
{
  return genulens::math::Interpolation::bilinear(nx, ny, F, xst, yst, dx, dy, xreq, yreq);
}
//---------------
void interp_xy_coeff(int nx, int ny, double *as, double xst, double yst, double dx, double dy, double xreq, double yreq) // just for this code
/* Return coefficients as[0], as[1], as[2], as[3] of bilinear interpolation.
 * Interpolated value is given by
 *   as[0]*F[ix][iy] + as[1]*F[ix+1][iy] + as[2]*F[ix][iy+1] + as[3]*F[ix+1][iy+1]  */
{
  genulens::math::Interpolation::bilinear_coefficients(nx, ny, as, xst, yst, dx, dy, xreq, yreq);
}

} // namespace genulens
