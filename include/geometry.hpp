#pragma once

#include <lib.hpp>

namespace hyperbolic {
    struct cosine {
        _float_t amp, phase, c; // amplitude, phase, and constant c
        cosine(_float_t amp, _float_t phase, _float_t c) : amp(amp), phase(phase), c(c) {}

        cosine operator+(const cosine &a) const;

        cosine operator-(const cosine &a) const;

        cosine operator*(_float_t a) const;

        _float_t operator()(_float_t x) const;

        [[nodiscard]] std::pair<_float_t, _float_t> zeros() const;
    };

    struct Bisector {
        _float_t numerator = 0;
        cosine denominator = {0, 0, 0};
        _float_t theta_start = 0, theta_end = 0;

        Bisector(pPoint s, pPoint t);

        void calc_phi();

        _float_t operator()(_float_t theta) const;

        [[nodiscard]] bool in_definition(_float_t theta) const;
    };
}