#include <cmath>
#include <utility>

#include <lib.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    _float_t clip(_float_t x);

    cosine cosine::operator +(const cosine& a) const {
        // the phase returned is always in [0, 2pi]
        _float_t x = cos(a.phase)*a.amp;
        x += cos(phase)*amp;
        _float_t y = sin(a.phase)*a.amp;
        y += sin(phase)*amp;

        _float_t result_amp = sqrt(pow(x, 2) + pow(y,2));
        _float_t result_phase = atan2(y, x);
        if (result_phase < 0) result_phase += 2*M_PI;
        return {result_amp, result_phase, a.c + c};
    }

    cosine cosine::operator -(const cosine& a) const{
        return *this + cosine(-a.amp, a.phase, -a.c);
    }

    cosine cosine::operator *(const _float_t a) const {
        return {a * amp, phase, a * c};
    }

    _float_t cosine::operator ()(_float_t x) const {
        return amp*cos(x + phase) + c;
    }

    [[nodiscard]] std::pair<_float_t, _float_t> cosine::zeros() const {
        // returns zeros of the cosine function a in the interval [0, 2 pi)
        _float_t z = acos(-c/amp);
        _float_t z1 = clip(z - phase);
        _float_t z2 = clip(2*M_PI - z - phase);
        return {z1, z2};
    }

    _float_t clip(_float_t x) {
        if (x < 0) return x + 2*M_PI;
        else if (x >= 2*M_PI) return x - 2*M_PI;
        return x;
    }

    _float_t calculate_beach_line_intersection(pSite siteS, pSite siteT, const _float_t r_sweep) {
        // returns angular coordinate of the intersection (s, t) under a sweep circle radius of r_sweep
        bool return_second = false;
        pPoint s = &siteS->point, t = &siteT->point;

        if (s->r < t->r) {
            std::swap(s, t);
            return_second = true;
        }

        cosine combined =
                cosine((cosh(r_sweep) - cosh(t->r))*sinh(s->r), 0, (cosh(t->r)-cosh(s->r))*sinh(r_sweep)) +
                cosine(-(cosh(r_sweep) - cosh(s->r))*sinh(t->r), s->theta-t->theta, 0);
        auto [z1, z2] = combined.zeros();
        if (z1 > z2) std::swap(z1, z2);

        if (return_second)
            return z2;
        return z1;
    }

    void Bisector::calc_phi() {
        // returns two angular coordinates between which the bisector is defined
        _float_t phi = acos(numerator/denominator.amp); // phi is in [0, pi/2]
        theta_start = clip(2*M_PI - phi - denominator.phase);
        theta_end = clip(phi - denominator.phase);
    }

    Bisector::Bisector(pPoint s, pPoint t) {
        if (s->r > t->r)
            std::swap(s, t);
        numerator = cosh(t->r) - cosh(s->r);
        denominator =
                cosine(sinh(t->r), -t->theta, 0) +
                cosine(-sinh(s->r), -s->theta, 0);
        calc_phi();
    }

    _float_t Bisector::operator()(_float_t theta) const {
        return atanh(numerator/denominator(theta));
    }

    [[nodiscard]] bool Bisector::in_definition(_float_t theta) const {
        theta = clip(theta - theta_start);
        _float_t end = clip(theta_end - theta_start);
        return theta <= end;
    }

    bool onActiveSide(rBeachLineElement a, rPoint p) {
        // checks if p is on the active side of the beach line intersection a

        rPoint s = a.first.point, t = a.second.point;
        _float_t outer_theta = (s.r >= t.r) ? s.theta : t.theta;
        _float_t p_theta = clip(p.theta - outer_theta);
        return (s.r >= t.r) ? p_theta <= M_PI : p_theta >= M_PI;
    }

    bool check_inside(Point& result, Bisector& rs, Bisector& st, _float_t z) {
        if (rs.in_definition(z) && st.in_definition(z)) {
            result.r = rs(z);
            result.theta = z;
            return true;
        }
        return false;
    }

    bool predict_circle_event(Point& result, pSite r, pSite s, pSite t) {
        if (!((r != s) && (s != t))) return false;

        // TODO: memorize already calculated values
        // TODO: use the hyperboloid model for more precise predictions

        Bisector rs(&r->point, &s->point);
        Bisector st(&s->point, &t->point);

        auto [z1, z2] = (rs.denominator*st.numerator - st.denominator*rs.numerator).zeros();

        if (check_inside(result, rs, st, z1)) return true;
        if (check_inside(result, rs, st, z2)) return true;
        return false;
    }

    bool predict_circle_event(Point& result, rBeachLineElement a, rBeachLineElement b) {
        if (!predict_circle_event(result, &a.first, &a.second, &b.second)) return false;
        return onActiveSide(a, result) && onActiveSide(b, result);
    }

    _float_t distance(rPoint s, rPoint t) {
        auto x = cosh(s.r)*cosh(t.r) - sinh(s.r)*cos(s.theta - t.theta)*sinh(t.r);
        return acosh(x);
    }
}