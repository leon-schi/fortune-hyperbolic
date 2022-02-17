#include <cmath>
#include <utility>

#include <lib.hpp>
#include <calculations.hpp>
#include <fortune.hpp>
#include <geometry.hpp>

#include <cassert>

/*
 * implementation of geometric calculations
 * */

namespace hyperbolic {

    _float_t calculate_beach_line_intersection(pSite siteS, pSite siteT, const double r_sweep) {
        // returns angular coordinate of the intersection (s, t) under a sweep circle radius of r_sweep
        bool return_second = false;
        pPoint s = &siteS->point, t = &siteT->point;

        if (s->r < t->r) {
            std::swap(s, t);
            return_second = true;
        }

        cosine<double> combined =
                cosine<double> ((cosh(r_sweep) - cosh(t->r)) * sinh(s->r), 0, (cosh(t->r) - cosh(s->r)) * sinh(r_sweep)) +
                cosine<double> (-(cosh(r_sweep) - cosh(s->r)) * sinh(t->r), s->theta - t->theta, 0);
        auto[z1, z2] = combined.zeros();
        if (z1 > z2) std::swap(z1, z2);

        z1 = clip(z1 + s->theta);
        z2 = clip(z2 + s->theta);

        return (return_second) ? z2 : z1;
    }

    bool assign_result_if_in_definiton(Point& result, Bisector<double>& rs, Bisector<double>& st, _float_t z) {
        if (rs.in_definition(z) && st.in_definition(z)) {
            result.r = rs(z);
            result.theta = z;
            return true;
        }
        return false;
    }

    bool assign_result_for_straight_bisector_if_exists(Point& result, Bisector<double>& straight, Bisector<double>& not_straight) {
        // assigns the coordinate of the intersection of straight and not_straight to result if it exists
        // straight is a straight line and not_straight is not

        assert(straight.is_straight);
        assert(!not_straight.is_straight);

        if (not_straight.in_definition(straight.straight_angle)) {
            result.r = not_straight(straight.straight_angle);
            result.theta = straight.straight_angle;
            return true;
        }
        return false;
    }

    bool predict_circle_event(Point& result, pSite r, pSite s, pSite t) {
        // TODO: maybe use hyperboloid model for predictions

        Bisector<double> rs(&r->point, &s->point);
        Bisector<double> st(&s->point, &t->point);

        // if one of the bisectors is a straight line, intersection calculation becomes straightforward
        if (rs.is_straight) {
            if (st.is_straight) {
                result = {0, 0};
                return true;
            }
            if (assign_result_for_straight_bisector_if_exists(result, rs, st)) return true;
        } else if (st.is_straight) {
            if (assign_result_for_straight_bisector_if_exists(result, st, rs)) return true;
        } else {
            auto [z1, z2] = (rs.denominator*st.numerator - st.denominator*rs.numerator).zeros();

            if (assign_result_if_in_definiton(result, rs, st, z1))
                return true;
            if (assign_result_if_in_definiton(result, rs, st, z2))
                return true;
        }
        return false;
    }
}