#pragma once

#include <cmath>
#include <utility>

#include <lib.hpp>
#include <datastructures.hpp>

/*
 * functions used for doing important geometric calculations
 * */

namespace hyperbolic {
    // transforms an angular coordinate into the range [0, 2*M_PI)
    _float_t clip(_float_t x);
    // calculates the hyperbolic distance between two points
    _float_t distance(rPoint s, rPoint t);
    // calculates the angular coordinate of the intersection of the beach curves of s and t at a sweep circle radius of r_sweep
    _float_t calculate_beach_line_intersection(pSite s, pSite  t, _float_t r_sweep);
}