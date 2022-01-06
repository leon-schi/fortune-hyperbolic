#pragma once

#include <cmath>
#include <utility>

#include <lib.hpp>
#include <geometry.hpp>
#include <datastructures.hpp>

namespace hyperbolic {
    _float_t distance(rPoint s, rPoint t);
    _float_t calculate_beach_line_intersection(pSite s, pSite  t, _float_t r_sweep);
    bool predict_circle_event(Point& result, rBeachLineElement a, rBeachLineElement b);
}