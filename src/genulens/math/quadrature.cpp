#include "genulens/math/quadrature.hpp"

namespace genulens::math {

int NewtonCotes::coefficients(int requested_order, double *locations, double *weights)
{
    int i = 0;
    int j = 0;
    int nmin = 1;
    const int order = (requested_order <= 1) ? 1 :
                      (requested_order <= 2) ? 2 :
                      (requested_order <= 4) ? 4 :
                      (requested_order <= 6) ? 6 :
                      (requested_order <= 8) ? 8 : 10;

    if (order == 1) {
        locations[i++] = 0.0;
        weights[j++] = 0.5;
        nmin = 1;
    }
    if (order == 2) {
        locations[i++] = 0.0, locations[i++] = 0.5, locations[i++] = 1.0;
        weights[j++] = 3.0 / 12, weights[j++] = 4.0 / 12, weights[j++] = 11.0 / 12;
        nmin = 3;
    }
    if (order == 4) {
        locations[i++] = 0.0,
        locations[i++] = 1.0 / 4, locations[i++] = 1.0 / 2, locations[i++] = 3.0 / 4,
        locations[i++] = 1.0, locations[i++] = 3.0 / 2, locations[i++] = 2.0,
        locations[i++] = 9.0 / 4, locations[i++] = 3.0;
        weights[j++] = 70. / 360,
        weights[j++] = 32. / 360, weights[j++] = 76. / 360, weights[j++] = 128. / 360,
        weights[j++] = 187. / 360, weights[j++] = 100. / 360,
        weights[j++] = 218. / 360, weights[j++] = 96. / 360, weights[j++] = 353. / 360;
        nmin = 7;
    }
    if (order == 6) {
        locations[i++] = 0.0,
        locations[i++] = 1.0 / 6, locations[i++] = 1.0 / 3, locations[i++] = 1.0 / 2,
        locations[i++] = 2.0 / 3, locations[i++] = 5.0 / 6, locations[i++] = 1.0,
        locations[i++] = 4.0 / 3, locations[i++] = 3.0 / 2, locations[i++] = 5.0 / 3,
        locations[i++] = 2.0, locations[i++] = 5.0 / 2, locations[i++] = 8.0 / 3,
        locations[i++] = 3.0, locations[i++] = 10.0 / 3, locations[i++] = 4.0,
        locations[i++] = 25.0 / 6, locations[i++] = 5.0;
        weights[j++] = 861. / 5040,
        weights[j++] = 216. / 5040, weights[j++] = 459. / 5040, weights[j++] = 920. / 5040,
        weights[j++] = 945. / 5040, weights[j++] = 1296. / 5040, weights[j++] = 2208. / 5040,
        weights[j++] = 162. / 5040, weights[j++] = 816. / 5040, weights[j++] = 567. / 5040,
        weights[j++] = 2955. / 5040, weights[j++] = 2008. / 5040, weights[j++] = 108. / 5040,
        weights[j++] = 3459. / 5040, weights[j++] = 999. / 5040, weights[j++] = 3662. / 5040,
        weights[j++] = 1080. / 5040, weights[j++] = 4999. / 5040;
        nmin = 11;
    }
    if (order == 8) {
        locations[i++] = 0.0,
        locations[i++] = 1.0 / 8, locations[i++] = 1.0 / 4, locations[i++] = 3.0 / 8,
        locations[i++] = 1.0 / 2, locations[i++] = 5.0 / 8, locations[i++] = 3.0 / 4,
        locations[i++] = 7.0 / 8, locations[i++] = 1.0, locations[i++] = 9.0 / 8,
        locations[i++] = 5.0 / 4, locations[i++] = 3.0 / 2, locations[i++] = 7.0 / 4,
        locations[i++] = 15.0 / 8, locations[i++] = 2.0, locations[i++] = 9.0 / 4,
        locations[i++] = 5.0 / 2, locations[i++] = 21.0 / 8, locations[i++] = 3.0,
        locations[i++] = 25.0 / 8, locations[i++] = 7.0 / 2, locations[i++] = 15.0 / 4,
        locations[i++] = 4.0, locations[i++] = 35.0 / 8, locations[i++] = 9.0 / 2,
        locations[i++] = 5.0, locations[i++] = 21.0 / 4, locations[i++] = 6.0,
        locations[i++] = 49.0 / 8, locations[i++] = 7.0;
        weights[j++] = 35604. / 226800,
        weights[j++] = 5888. / 226800, weights[j++] = 10848. / 226800, weights[j++] = 28160. / 226800,
        weights[j++] = 17156. / 226800, weights[j++] = 39936. / 226800, weights[j++] = 52608. / 226800,
        weights[j++] = 47104. / 226800, weights[j++] = 43213. / 226800, weights[j++] = 31488. / 226800,
        weights[j++] = 16352. / 226800, weights[j++] = 20940. / 226800, weights[j++] = 5280. / 226800,
        weights[j++] = 83968. / 226800, weights[j++] = 31410. / 226800, weights[j++] = 60192. / 226800,
        weights[j++] = 19284. / 226800, weights[j++] = 91136. / 226800, weights[j++] = 103575. / 226800,
        weights[j++] = 52480. / 226800, weights[j++] = -8228. / 226800, weights[j++] = 58336. / 226800,
        weights[j++] = 99196. / 226800, weights[j++] = 102912. / 226800, weights[j++] = -5568. / 226800,
        weights[j++] = 184153. / 226800, weights[j++] = 28832. / 226800, weights[j++] = 177718. / 226800,
        weights[j++] = 41216. / 226800, weights[j++] = 225811. / 226800;
        nmin = 15;
    }
    if (order == 10) {
        locations[i++] = 0.0,
        locations[i++] = 1.0 / 10, locations[i++] = 1.0 / 5, locations[i++] = 3.0 / 10,
        locations[i++] = 2.0 / 5, locations[i++] = 1.0 / 2, locations[i++] = 3.0 / 5,
        locations[i++] = 7.0 / 10, locations[i++] = 4.0 / 5, locations[i++] = 9.0 / 10,
        locations[i++] = 1.0, locations[i++] = 6.0 / 5, locations[i++] = 7.0 / 5,
        locations[i++] = 3.0 / 2, locations[i++] = 8.0 / 5, locations[i++] = 9.0 / 5,
        locations[i++] = 2.0, locations[i++] = 21.0 / 10, locations[i++] = 12.0 / 5,
        locations[i++] = 5.0 / 2, locations[i++] = 27.0 / 10, locations[i++] = 14.0 / 5,
        locations[i++] = 3.0, locations[i++] = 16.0 / 5, locations[i++] = 7.0 / 2,
        locations[i++] = 18.0 / 5, locations[i++] = 4.0, locations[i++] = 21.0 / 5,
        locations[i++] = 9.0 / 2, locations[i++] = 24.0 / 5, locations[i++] = 49.0 / 10,
        locations[i++] = 5.0, locations[i++] = 27.0 / 5, locations[i++] = 28.0 / 5,
        locations[i++] = 6.0, locations[i++] = 63.0 / 10, locations[i++] = 32.0 / 5,
        locations[i++] = 7.0, locations[i++] = 36.0 / 5, locations[i++] = 8.0,
        locations[i++] = 81.0 / 10, locations[i++] = 9.0;
        weights[j++] = 883685. / 5987520,
        weights[j++] = 106300. / 5987520, weights[j++] = 164075. / 5987520, weights[j++] = 591300. / 5987520,
        weights[j++] = 67600. / 5987520, weights[j++] = 958868. / 5987520, weights[j++] = 776475. / 5987520,
        weights[j++] = 1016500. / 5987520, weights[j++] = 86675. / 5987520, weights[j++] = 1880200. / 5987520,
        weights[j++] = 1851848. / 5987520, weights[j++] = -504300. / 5987520, weights[j++] = 205125. / 5987520,
        weights[j++] = 2644104. / 5987520, weights[j++] = -1527450. / 5987520, weights[j++] = 628625. / 5987520,
        weights[j++] = 1177276. / 5987520, weights[j++] = 2724000. / 5987520, weights[j++] = -571875. / 5987520,
        weights[j++] = 2136840. / 5987520, weights[j++] = 2770500. / 5987520, weights[j++] = -734250. / 5987520,
        weights[j++] = 4772079. / 5987520, weights[j++] = -2278500. / 5987520, weights[j++] = 4353576. / 5987520,
        weights[j++] = -3483050. / 5987520, weights[j++] = 4097507. / 5987520, weights[j++] = -189450. / 5987520,
        weights[j++] = 4377812. / 5987520, weights[j++] = -2375550. / 5987520, weights[j++] = 1906800. / 5987520,
        weights[j++] = 5210935. / 5987520, weights[j++] = -1707150. / 5987520, weights[j++] = 1839525. / 5987520,
        weights[j++] = 2621502. / 5987520, weights[j++] = 3195700. / 5987520, weights[j++] = -388200. / 5987520,
        weights[j++] = 5361569. / 5987520, weights[j++] = 413675. / 5987520, weights[j++] = 4892386. / 5987520,
        weights[j++] = 956700. / 5987520, weights[j++] = 5971453. / 5987520;
        nmin = 19;
    }
    return nmin;
}

} // namespace genulens::math
