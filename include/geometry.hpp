#pragma once

#include <cmath>
#include <vector>

using std::vector, std::max, std::min;

namespace hyperbolic {
    typedef double _float_t;

    template<typename _float_T>  struct _Point;
    template<typename _float_T> struct HyperboloidVec;

    /*
     * Most basic class for representing a point in polar coordinates
     * */
    template<typename _float_T>
    struct _Point {
        _float_T r = 0, theta = 0;
        _Point() = default;
        _Point(_float_T r, _float_T theta) : r(r), theta(theta) {}
        explicit _Point(const HyperboloidVec<_float_T>& p);
        bool operator < (const _Point<_float_T>& p) const {
            return r < p.r;
        }
    };
    typedef _Point<double> Point;
    typedef const Point* pPoint;
    typedef const Point& rPoint;

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
}