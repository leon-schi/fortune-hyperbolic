#pragma once

#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <type_traits>
#include <unordered_map>

#include <fortune-hyperbolic/geometry.hpp>

using std::is_base_of, std::unique_ptr, std::shared_ptr, std::vector, std::priority_queue, std::size_t, std::string, std::to_string, std::hash, std::unordered_map;

namespace hyperbolic {
    template<typename _float_T>
    class BeachLineElement;

    /**
     * base class for the events handled in the algorithm
     * */
    template<typename _float_T>
    class Event {
    public:
        const _float_T r {0};
        Event(_float_T r) : r(r) {};
        virtual ~Event() = default;
    };

    template<typename _float_T>
    class SiteEvent : public Event<_float_T> {
    public:
        Site<_float_T>& site;
        explicit SiteEvent(Site<_float_T>& s) : Event<_float_T>(s.point.r), site(s) {};
    };

    template<typename _float_T>
    class CircleEvent : public Event<_float_T> {
    private:
        bool valid = true;
    public:
        // references to the beach line elements that cause the circle event
        BeachLineElement<_float_T> &first, &second;
        // the point at which a Voronoi vertex is created when the event happens
        const Point<_float_T> center;
        CircleEvent(BeachLineElement<_float_T>& first, BeachLineElement<_float_T>& second, Point<_float_T> center, _float_T r)
                : Event<_float_T>(r + center.r), first(first), second(second), center(center) {};

        void invalidate() {
            valid = false;
        }
        [[nodiscard]] bool isValid() const {
            return valid;
        }
    };

    // comparison of events based on their radius
    template <class T, typename _float_T>
    struct CmpEvent {
        static_assert(is_base_of<Event<_float_T>, T>::value, "template type T must inherit from Event");
        static bool greater_equal(const Event<_float_T>* e1, const Event<_float_T>* e2) {
            return e1->r > e2->r;
        }

        bool operator()(const shared_ptr<T>& e1, const shared_ptr<T> &e2) {
            return greater_equal(e1.get(), e2.get());
        }
    };

    // queue of events
    template <class T, typename _float_T>
    class EventQueue : public priority_queue<shared_ptr<T>, vector<shared_ptr<T>>, CmpEvent<T, _float_T>> {
    public:
        T* getTop() {
            return (this->empty()) ? nullptr : this->top().get();
        }
    };

    using ull = unsigned long long;
    // struct identifying a set of three sites based on their IDs. Used for caching the prediction of circle events
    struct SiteTriple {
        ull ID1, ID2, ID3;
        SiteTriple(ull a, ull b, ull c) {
            vector<unsigned long long> ids = {a, b, c};
            std::sort(ids.begin(), ids.end());
            ID1 = ids[0], ID2 = ids[1], ID3 = ids[2];
        }
    };
    struct CmpSiteTriple {
        bool operator ()(const SiteTriple& a, const SiteTriple& b) const {
            return a.ID1 == b.ID1 && a.ID2 == b.ID2 && a.ID3 == b.ID3;
        }
    };
    struct HashSiteTriple {
        size_t operator ()(const SiteTriple& a) const {
            string s = to_string(a.ID1) + to_string(a.ID2) + to_string(a.ID3);
            return hash<string>{}(s);
        }
    };

    template<typename _float_T>
    using SiteTripleMap = unordered_map<SiteTriple, Point<_float_T>, HashSiteTriple, CmpSiteTriple>;


}