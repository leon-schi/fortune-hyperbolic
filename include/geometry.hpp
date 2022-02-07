#pragma once

#include <cmath>
#include <vector>

#include <lib.hpp>

using std::vector, std::max, std::min;

namespace hyperbolic {
    typedef double _float_t;

    struct Point;
    struct HyperboloidVec;

    /*
     * Most basic class for representing a point in polar coordinates
     * */
    struct Point {
        _float_t r = 0, theta = 0;
        Point() = default;
        Point(_float_t r, _float_t theta) : r(r), theta(theta) {}
        explicit Point(const HyperboloidVec& p);
        bool operator < (const Point& p) const {
            return r < p.r;
        }
    };
    typedef const Point* pPoint;
    typedef const Point& rPoint;

    // represents a vector in the hyperboloid model
    struct HyperboloidVec {
        _float_t x, y, z;
        HyperboloidVec() = default;
        explicit HyperboloidVec(const Point& p);
        HyperboloidVec(_float_t x, _float_t y, _float_t z) : x(x), y(y), z(z) {};
        HyperboloidVec operator+(const HyperboloidVec& v) const;
        HyperboloidVec operator-(const HyperboloidVec& v) const;
        HyperboloidVec operator*(double a) const;
        _float_t dot(const HyperboloidVec& v) const;
        HyperboloidVec cross(HyperboloidVec& v) const;
        void normalize();
    };

    // represents a Bisector as a parameterized line in the hyperboloid model
    struct HyperboloidBisector {
        HyperboloidVec u, v;
        void find_u(pPoint a, pPoint b);
        HyperboloidBisector(pPoint a, pPoint b);
    };

    /*
     * Sites for internal use that assign an ID to a point which is used for caching
     * */
    struct Site {
        rPoint point;
        const unsigned long long ID;
        Site(rPoint point, unsigned long long id) : point(point), ID(id) {};
    };
    typedef const Site & rSite;
    typedef const Site * pSite;

    /*
     * struct representing a cosine function for geometric calculations
     * */
    struct cosine {
        _float_t amp, phase, c; // amplitude, phase, and constant c
        cosine(_float_t amp, _float_t phase, _float_t c) : amp(amp), phase(phase), c(c) {}

        cosine operator+(const cosine &a) const;

        cosine operator-(const cosine &a) const;

        cosine operator*(_float_t a) const;

        _float_t operator()(_float_t x) const;

        [[nodiscard]] std::pair<_float_t, _float_t> zeros() const;
    };

    /*
     * struct representing a line bisector between two points used for geometric calculations and drawing
     * */
    struct Bisector {
        /*
         * the function is represented as r(theta) = atanh(numerator / denominator(theta))
         * */
        _float_t numerator = 0;
        cosine denominator = {0, 0, 0};

        // marks whether the object represents the special case of a straight line and if so at which angular coordinate it is placed
        bool is_straight = false;
        _float_t straight_angle = 0;


        // the angular coordinates between which the bisector is defined in case it is not a straight bisector
        _float_t theta_start = 0, theta_end = 0;
        void calc_definition();


        // returns the radial coordinate at angular coordinate theta
        _float_t operator()(_float_t theta) const;

        // constructs the object from the two points that define it
        Bisector(pPoint s, pPoint t);

        // checks whether the bisector is defined at the angular coordinate theta
        [[nodiscard]] bool in_definition(_float_t theta) const;
    };
}