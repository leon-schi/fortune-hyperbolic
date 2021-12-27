#include <cmath>
#include <utility>

#include <lib.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    struct cosine {
        _float_t amp, phase, c; // amplitude, phase in the interval [-pi, pi], and constant c
    };

    cosine add_cosines(cosine a, cosine b) {
        _float_t x = cos(a.phase)*a.amp;
        x += cos(b.phase)*b.amp;
        _float_t y = sin(a.phase)*a.amp;
        y += sin(b.phase)*b.amp;

        _float_t amp = sqrt(pow(x, 2) + pow(y,2));
        _float_t phase = atan2(y, x);
        return {amp, phase, a.c + b.c};
    }

    _float_t clip(_float_t x) {
        if (x < 0) return x + 2*M_PI;
        else if (x >= 2*M_PI) return x - 2*M_PI;
        return x;
    }

    std::pair<_float_t, _float_t> zeros(cosine& a) {
        // returns zeros of the cosine function a in the interval [0, 2 pi) ordered increasingly
        _float_t z = acos(-a.c/a.amp);
        _float_t z1 = clip(z - a.phase);
        _float_t z2 = clip(2*M_PI - z - a.phase);

        if (z1 > z2)
            std::swap(z1, z2);

        return {z1, z2};
    }

    _float_t calculate_beach_line_intersection(pSite s, pSite t, const _float_t r_sweep) {
        // returns angular coordinate of the intersection (s, t) under a sweep circle radius of r_sweep
        bool return_second = false;
        if (s->r < t->r) {
            std::swap(s, t);
            return_second = true;
        }

        cosine combined = add_cosines(
                {cosh(s->r)*sinh(t->r), s->theta-t->theta, -cosh(s->r)*sinh(r_sweep)},
                {-cosh(t->r)*sinh(s->r), 0, cosh(t->r)*sinh(r_sweep)});
        auto [z1, z2] = zeros(combined);

        if (return_second)
            return z2;
        return z1;
    }
}