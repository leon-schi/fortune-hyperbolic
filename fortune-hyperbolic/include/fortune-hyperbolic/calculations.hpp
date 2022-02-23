#pragma once

#include <cmath>
#include <utility>

#include <fortune-hyperbolic/geometry.hpp>
#include <fortune-hyperbolic/datastructures.hpp>

/*
 * functions used for doing important geometric calculations
 * */

namespace hyperbolic {

    // represents a vector in the hyperboloid model
    template<typename _float_T>
    struct HyperboloidVec {
        _float_T x, y, z;
        HyperboloidVec() = default;
        explicit HyperboloidVec(const Point<_float_T>& p) {
            x = sinh(p.r) * cos(p.theta);
            y = sinh(p.r) * sin(p.theta);
            z = cosh(p.r);
        };
        HyperboloidVec(_float_T x, _float_T y, _float_T z) : x(x), y(y), z(z) {};
        HyperboloidVec operator+(const HyperboloidVec& v) const {
            return HyperboloidVec(x + v.x, y + v.y, z + v.z);
        };
        HyperboloidVec operator*(double a) const {
            return HyperboloidVec(a*x, a*y, a*z);
        };
        HyperboloidVec operator-(const HyperboloidVec& v) const {
            return *this + v*(-1);
        };
        _float_T dot(const HyperboloidVec& v) const {
            return x*v.x + y*v.y - z*v.z;
        };
        HyperboloidVec cross(HyperboloidVec& v) const {
            return HyperboloidVec(y*v.z - v.y*z, v.x*z - x*v.z, y*v.x - v.y*x);
        };
        void normalize() {
            _float_T a_sq = z*z - x*x - y*y;
            if (a_sq < 0) a_sq *= -1;
            _float_T a = sqrt(a_sq);
            if (z < 0) a = -a;
            x /= a; y /= a; z/= a;
        };
    };

    // represents a Bisector as a parameterized line in the hyperboloid model
    template<typename _float_T>
    struct HyperboloidBisector {
        HyperboloidVec<_float_T> u, v;
        void find_u(Point<_float_T>* a, Point<_float_T>* b) {
            if (a->r > b->r) std::swap(a, b);
            HyperboloidVec<_float_T> v_a(*a), v_b(*b);

            HyperboloidVec<_float_T> p_b(-v_b.y, v_b.x, 0);
            HyperboloidVec<_float_T> p_ab = v_a - v_b;

            u = p_b.cross(p_ab);
            u.normalize();
        };
        HyperboloidBisector(Point<_float_T>& a, Point<_float_T>& b) {
            find_u(&a, &b);
            HyperboloidVec<_float_T> v_a(a), v_b(b);
            HyperboloidVec<_float_T> v_ab = v_b - v_a;
            v = u.cross(v_ab);
            v.normalize();
        };
    };

    // transforms an angular coordinate into the range [0, 2*M_PI)
    template <typename _float_T>
    _float_T clip(_float_T x) {
        if (x < 0) return x + 2*M_PI;
        else if (x >= 2*M_PI) return x - 2*M_PI;
        return x;
    }

    // calculates the hyperbolic distance between two points
    template <typename _float_T>
    _float_T distance(const Point<_float_T>& s, const Point<_float_T>& t) {
        auto x = cosh(s.r)*cosh(t.r) - sinh(s.r)*cos(s.theta - t.theta)*sinh(t.r);
        return acosh(x);
    };

    /*
     * struct representing a cosine function for geometric calculations
     * */
    template <typename _float_T>
    struct cosine {
        _float_T amp, phase, c; // amplitude, phase, and constant c
        cosine(_float_T amp, _float_T phase, _float_T c) : amp(amp), phase(phase), c(c) {}

        cosine<_float_T> operator+(const cosine<_float_T>& a) const {
            // the phase returned is always in [0, 2pi]
            _float_T x = cos(a.phase)*a.amp;
            x += cos(phase)*amp;
            _float_T y = sin(a.phase)*a.amp;
            y += sin(phase)*amp;

            _float_T result_amp = sqrt(pow(x, 2) + pow(y,2));
            _float_T result_phase = atan2(y, x);
            if (result_phase < 0) result_phase += 2*M_PI;
            return {result_amp, result_phase, a.c + c};
        }

        cosine<_float_T> operator-(const cosine<_float_T>& a) const{
            return *this + cosine(-a.amp, a.phase, -a.c);
        }

        cosine<_float_T> operator *(const _float_T a) const {
            return {a * amp, phase, a * c};
        }

        _float_T operator ()(_float_T x) const {
            return amp*cos(x + phase) + c;
        }

        [[nodiscard]] std::pair<_float_T, _float_T> zeros() const {
            // returns zeros of the cosine function a in the interval [0, 2 pi)
            _float_T z = acos(-c/amp);
            _float_T z1 = clip<_float_T>(z - phase);
            _float_T z2 = clip<_float_T>(2*M_PI - z - phase);
            return {z1, z2};
        }
    };

    /**
     * struct representing a line bisector between two points used for geometric calculations and drawing
     * the function is represented as r(theta) = atanh(numerator / denominator(theta))
     * */
    template <typename _float_T>
    struct Bisector {
        _float_T numerator = 0;
        cosine<_float_T> denominator = {0, 0, 0};

        // the angular coordinates between which the bisector is defined in case it is not a straight bisector
        _float_T theta_start = 0, theta_end = 0;

        // marks whether the object represents the special case of a straight line and if so at which angular coordinate it is placed
        bool is_straight = false;
        _float_T straight_angle = 0;

        void calc_definition() {
            // returns two angular coordinates between which the bisector is defined
            _float_T phi = acos(numerator/denominator.amp); // phi is in [0, pi/2]
            theta_start = clip<_float_T>(2*M_PI - phi - denominator.phase);
            theta_end = clip<_float_T>(phi - denominator.phase);
        }

        // constructs the object from the two points that define it
        Bisector(const Point<_float_T>* s, const Point<_float_T>* t) {
            if (s->r > t->r)
                std::swap(s, t);
            numerator = cosh(t->r) - cosh(s->r);

            if (numerator == 0.0) {
                is_straight = true;
                straight_angle = (s->theta + t->theta) / 2.0;
            } else {
                is_straight = false;
                denominator =
                        cosine<_float_T>(sinh(t->r), -t->theta, 0) +
                        cosine<_float_T>(-sinh(s->r), -s->theta, 0);
                calc_definition();
            }
        }

        // returns the radial coordinate at angular coordinate theta
        _float_T operator()(_float_T theta) const {
            return atanh(numerator/denominator(theta));
        }

        // checks whether the bisector is defined at the angular coordinate theta
        [[nodiscard]] bool in_definition(_float_T theta) const {
            theta = clip<_float_T>(theta - theta_start);
            _float_T end = clip<_float_T>(theta_end - theta_start);
            return theta <= end;
        }
    };

    template<typename _float_T>
    Point<_float_T>::Point(const HyperboloidVec<_float_T>& p) {
        r = acosh(p.z); theta = clip(atan2(p.y, p.x));
    }
}