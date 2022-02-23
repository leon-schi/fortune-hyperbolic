#pragma once

#include <fortune-hyperbolic/geometry.hpp>
#include <fortune-hyperbolic/calculations.hpp>
#include <fortune-hyperbolic/datastructures.hpp>

namespace hyperbolic {
    template<typename _float_T>
    class Kernel {
    public:
        /*
         * checks whether the angle theta is between the reference angle and the Beach Line intersection defined by the tuple of sites (s,t)
         * */
        virtual bool before (
                _float_T theta, _float_T reference_angle,
                Point<_float_T>& p_s, Point<_float_T>& p_t, _float_T r_sweep) = 0;

        /*
         * predict Circle Event
         * */
        virtual bool predict_circle_event(
                Point<_float_T>& result,
                Site<_float_T>& r, Site<_float_T>& s, Site<_float_T>& t)  = 0;
    };

    template<typename _float_T>
    class FullNativeKernel: public Kernel<_float_T> {
        SiteTripleMap<_float_T> circleEventCache;
    public:
        /*
         * returns true if theta is between reference_angle and the beach line intersection (s,t) in ccw direction
         * */
        bool before (
                _float_T theta, _float_T reference_angle,
                Point<_float_T>& p_s, Point<_float_T>& p_t, _float_T r_sweep) {
            Point<_float_T> *s = &p_s, *t = &p_t;

            auto r_sweep_internal = static_cast<_float_T>(r_sweep);

            bool use_second = false;
            if (s->r < t->r) {
                std::swap(s, t);
                use_second = true;
            }

            cosine<_float_T> combined =
                    cosine<_float_T> ((cosh(r_sweep_internal) - cosh(t->r)) * sinh(s->r), 0, (cosh(t->r) - cosh(s->r)) * sinh(r_sweep_internal)) +
                    cosine<_float_T> (-(cosh(r_sweep_internal) - cosh(s->r)) * sinh(t->r), s->theta - t->theta, 0);
            auto[z1, z2] = combined.zeros();
            if (z1 > z2) std::swap(z1, z2);

            z1 = clip<_float_T>(z1 + s->theta);
            z2 = clip<_float_T>(z2 + s->theta);
            auto z = (use_second) ? z2 : z1;

            z = clip<_float_T>(z - reference_angle);
            return (z > theta);
        };

        bool assign_result_if_in_definiton(Point<_float_T>& result, Bisector<_float_T>& rs, Bisector<_float_T>& st, _float_T z) {
            if (rs.in_definition(z) && st.in_definition(z)) {
                result.r = rs(z);
                result.theta = z;
                return true;
            }
            return false;
        }

        bool assign_result_for_straight_bisector_if_exists(Point<_float_T>& result, Bisector<_float_T>& straight, Bisector<_float_T>& not_straight) {
            // assigns the coordinate of the intersection of straight and not_straight to result if it exists
            // straight is a straight line and not_straight is not

            if (not_straight.in_definition(straight.straight_angle)) {
                result.r = not_straight(straight.straight_angle);
                result.theta = straight.straight_angle;
                return true;
            }
            return false;
        }

        bool calculate_circle_event_center(Point<_float_T>& result, const Point<_float_T>& r, const Point<_float_T>& s, const Point<_float_T>& t) {
            Bisector<_float_T> rs(&r, &s);
            Bisector<_float_T> st(&s, &t);

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
                auto zeros = (rs.denominator*st.numerator - st.denominator*rs.numerator).zeros();

                if (assign_result_if_in_definiton(result, rs, st, zeros.first))
                    return true;
                if (assign_result_if_in_definiton(result, rs, st, zeros.second))
                    return true;
            }
            return false;
        }

        // checks if p is on the active side of the beach line intersection a
        bool on_active_site(Point<_float_T>& s, Point<_float_T>& t, Point<_float_T>& p) {
            if (p.r == 0.0) return true;
            _float_t outer_theta = (s.r >= t.r) ? s.theta : t.theta;
            _float_t p_theta = clip<_float_T>(p.theta - outer_theta);
            return (s.r >= t.r) ? p_theta <= M_PI : p_theta >= M_PI;
        }

        bool predict_circle_event(
                Point<_float_T>& result,
                Site<_float_T>& r, Site<_float_T>& s, Site<_float_T>& t) {

            auto siteTriple = SiteTriple(r.ID, s.ID, t.ID);
            if (siteTriple.ID1 == siteTriple.ID2 || siteTriple.ID2 == siteTriple.ID3) return false;

            auto it = circleEventCache.find(siteTriple);
            if (it != circleEventCache.end()) {
                // use cached value
                result = it->second;
            } else {
                // calculate point and cache it if existent
                if (calculate_circle_event_center(result, r.point, s.point, t.point))
                    circleEventCache[siteTriple] = result;
                else return false;
            }

            return (on_active_site(r.point, s.point, result) &&
                    on_active_site(s.point, t.point, result));
        };
    };
}