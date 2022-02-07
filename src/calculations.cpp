#include <cmath>
#include <utility>

#include <lib.hpp>
#include <calculations.hpp>
#include <fortune.hpp>
#include <geometry.hpp>

#include <cassert>

/*
 * implementation of geometric calculations
 * */

namespace hyperbolic {
    HyperboloidVec HyperboloidVec::cross(HyperboloidVec& v) const {
        return HyperboloidVec(y*v.z - v.y*z, v.x*z - x*v.z, y*v.x - v.y*x);
    }

    void HyperboloidVec::normalize() {
        _float_t a_sq = z*z - x*x - y*y;
        if (a_sq < 0) a_sq *= -1;
        _float_t a = sqrt(a_sq);
        if (z < 0) a = -a;
        x /= a; y /= a; z/= a;
    }

    HyperboloidVec::HyperboloidVec(const Point& p) {
        x = sinh(p.r) * cos(p.theta);
        y = sinh(p.r) * sin(p.theta);
        z = cosh(p.r);
    }

    HyperboloidVec HyperboloidVec::operator+(const HyperboloidVec& v) const {
        return HyperboloidVec(x + v.x, y + v.y, z + v.z);
    }

    HyperboloidVec HyperboloidVec::operator-(const HyperboloidVec& v) const {
        return HyperboloidVec(x - v.x, y - v.y, z - v.z);
    }

    HyperboloidVec HyperboloidVec::operator*(double a) const {
        return HyperboloidVec(a*x, a*y, a*z);
    }

    _float_t HyperboloidVec::dot(const HyperboloidVec& v) const {
        return x*v.x + y*v.y - z*v.z;
    }

    void HyperboloidBisector::find_u(pPoint a, pPoint b) {
        if (a->r > b->r) std::swap(a, b);
        HyperboloidVec v_a(*a), v_b(*b);

        HyperboloidVec p_b(-v_b.y, v_b.x, 0);
        HyperboloidVec p_ab = v_a - v_b;

        u = p_b.cross(p_ab);
        u.normalize();
    }

    HyperboloidBisector::HyperboloidBisector(pPoint a, pPoint b) {
        find_u(a,b);
        HyperboloidVec v_a(*a), v_b(*b);
        HyperboloidVec v_ab = v_b - v_a;
        v = u.cross(v_ab);
        v.normalize();
    }

    Point::Point(const HyperboloidVec& p) {
        r = acosh(p.z); theta = clip(atan2(p.y, p.x));
    }

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
                cosine((cosh(r_sweep) - cosh(t->r)) * sinh(s->r), 0, (cosh(t->r) - cosh(s->r)) * sinh(r_sweep)) +
                cosine(-(cosh(r_sweep) - cosh(s->r)) * sinh(t->r), s->theta - t->theta, 0);
        auto[z1, z2] = combined.zeros();
        if (z1 > z2) std::swap(z1, z2);

        z1 = clip(z1 + s->theta);
        z2 = clip(z2 + s->theta);

        return (return_second) ? z2 : z1;
    }

    void Bisector::calc_definition() {
        // returns two angular coordinates between which the bisector is defined
        _float_t phi = acos(numerator/denominator.amp); // phi is in [0, pi/2]
        theta_start = clip(2*M_PI - phi - denominator.phase);
        theta_end = clip(phi - denominator.phase);
    }

    Bisector::Bisector(pPoint s, pPoint t) {
        if (s->r > t->r)
            std::swap(s, t);
        numerator = cosh(t->r) - cosh(s->r);

        if (numerator == 0.0) {
            is_straight = true;
            straight_angle = (s->theta + t->theta) / 2.0;
        } else {
            is_straight = false;
            denominator =
                    cosine(sinh(t->r), -t->theta, 0) +
                    cosine(-sinh(s->r), -s->theta, 0);
            calc_definition();
        }
    }

    _float_t Bisector::operator()(_float_t theta) const {
        return atanh(numerator/denominator(theta));
    }

    [[nodiscard]] bool Bisector::in_definition(_float_t theta) const {
        theta = clip(theta - theta_start);
        _float_t end = clip(theta_end - theta_start);
        return theta <= end;
    }

    // checks if p is on the active side of the beach line intersection a
    bool onActiveSide(rBeachLineElement a, rPoint p) {
        if (p.r == 0.0) return true;

        rPoint s = a.first.point, t = a.second.point;
        _float_t outer_theta = (s.r >= t.r) ? s.theta : t.theta;
        _float_t p_theta = clip(p.theta - outer_theta);
        return (s.r >= t.r) ? p_theta <= M_PI : p_theta >= M_PI;
    }

    _float_t distance(rPoint s, rPoint t) {
        auto x = cosh(s.r)*cosh(t.r) - sinh(s.r)*cos(s.theta - t.theta)*sinh(t.r);
        return acosh(x);
    }

    bool assign_result_if_in_definiton(Point& result, Bisector& rs, Bisector& st, _float_t z) {
        if (rs.in_definition(z) && st.in_definition(z)) {
            result.r = rs(z);
            result.theta = z;
            return true;
        }
        return false;
    }

    bool assign_result_for_straight_bisector_if_exists(Point& result, Bisector& straight, Bisector& not_straight) {
        // assigns the coordinate of the intersection of straight and not_straight to result if it exists
        // straight is a straight line and not_straight is not

        assert(straight.is_straight);
        assert(!not_straight.is_straight);

        if (not_straight.in_definition(straight.straight_angle)) {
            result.r = not_straight(straight.straight_angle);
            result.theta = straight.straight_angle;
            return true;
        }
        return false;
    }

    bool predict_circle_event(Point& result, pSite r, pSite s, pSite t) {
        // TODO: maybe use hyperboloid model for predictions

        Bisector rs(&r->point, &s->point);
        Bisector st(&s->point, &t->point);

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
            auto [z1, z2] = (rs.denominator*st.numerator - st.denominator*rs.numerator).zeros();

            if (assign_result_if_in_definiton(result, rs, st, z1))
                return true;
            if (assign_result_if_in_definiton(result, rs, st, z2))
                return true;
        }
        return false;
    }

    bool FortuneHyperbolicImplementation::predictCircleEvent(Point& result, rBeachLineElement a, rBeachLineElement b) {
        rSite r = a.first, s = a.second, t = b.second;
        auto siteTriple = SiteTriple(r, s, t);
        if (siteTriple.ID1 == siteTriple.ID2 || siteTriple.ID2 == siteTriple.ID3) return false;

        auto it = circleEventCache.find(siteTriple);
        if (it != circleEventCache.end()) {
            // use cached value
            result = it->second;
        } else {
            // calculate point and cache it if existent
            if (hyperbolic::predict_circle_event(result, &r, &s, &t))
                circleEventCache[siteTriple] = result;
            else return false;
        }

        return onActiveSide(a, result) && onActiveSide(b, result);
    }
}