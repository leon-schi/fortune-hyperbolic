#pragma once

#include <geometry.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    template<typename _float_T_internal, typename _float_T_external>
    class Kernel {
    public:
        /*
         * checks whether the angle theta is between the reference angle and the Beach Line intersection defined by the tuple of sites (a,b)
         * */
        virtual bool before(
                _float_T_external theta, _float_T_external reference_angle,
                const _Point<_float_T_external>& SiteS, const _Point<_float_T_external>& SiteT,
                _float_T_external r_sweep) = 0;

        /*
         * predict Circle Event
         * */
        virtual bool predict_circle_event(
                _Point<_float_T_external>& result,
                const _Point<_float_T_external>& a, const _Point<_float_T_external>& b, const _Point<_float_T_external>& c) = 0;
    };

    template<typename _float_T_internal, typename _float_T_external>
    class FullNativeKernel: public Kernel<_float_T_internal, _float_T_external> {
        typedef _Point<_float_T_external> Point_ext;
        typedef _Point<_float_T_internal> Point_int;

        static Point_int convert(Point_ext p) {
            return Point_int(static_cast<_float_T_internal>(p.r), static_cast<_float_T_internal>(p.theta));
        }

        static long double sinh(long double x) {
            return sinhl(x);
        }

        static long double cosh(long double x) {
            return coshl(x);
        }

    public:
        bool before(
                _float_T_external theta, _float_T_external reference_angle,
                const _Point<_float_T_external>& SiteS, const _Point<_float_T_external>& SiteT,
                _float_T_external r_sweep) override {
            // find the angular coordinate of the intersection (s, t) under a sweep circle radius of r_sweep
            Point_int p_s = convert(SiteS), p_t = convert(SiteT);
            Point_int *s = &p_s, *t = &p_t;

            auto r_sweep_internal = static_cast<_float_T_internal>(r_sweep);

            bool use_second = false;
            if (s->r < t->r) {
                std::swap(s, t);
                use_second = true;
            }

            cosine<_float_T_internal> combined =
                    cosine<_float_T_internal> ((cosh(r_sweep_internal) - cosh(t->r)) * sinh(s->r), 0, (cosh(t->r) - cosh(s->r)) * sinh(r_sweep_internal)) +
                    cosine<_float_T_internal> (-(cosh(r_sweep_internal) - cosh(s->r)) * sinh(t->r), s->theta - t->theta, 0);
            auto[z1, z2] = combined.zeros();
            if (z1 > z2) std::swap(z1, z2);

            z1 = clip<_float_T_internal>(z1 + s->theta);
            z2 = clip<_float_T_internal>(z2 + s->theta);
            auto z = (use_second) ? z2 : z1;

            z = clip<_float_T_internal>(z - reference_angle);
            return (z <= theta);
        };

        bool predict_circle_event(
                _Point<_float_T_external>& result,
                const _Point<_float_T_external>& a, const _Point<_float_T_external>& b, const _Point<_float_T_external>& c) override {
            return false;
        };
    };
}