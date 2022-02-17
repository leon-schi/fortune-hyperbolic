#pragma once

#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <type_traits>
#include <unordered_map>

#include <lib.hpp>

using std::is_base_of, std::unique_ptr, std::shared_ptr, std::vector, std::priority_queue, std::size_t, std::string, std::to_string, std::hash, std::unordered_map;

namespace hyperbolic {
    class BeachLineElement;
    typedef BeachLineElement* pBeachLineElement;
    typedef BeachLineElement& rBeachLineElement;

    /*
     * base class for the events handled in the algorithm
     * */
    class Event {
    public:
        const _float_t r {0};
        explicit Event(_float_t r) : r(r) {};
        virtual ~Event() = default;
    };

    class SiteEvent : public Event {
    public:
        rSite site;
        explicit SiteEvent(const Site& s) : Event(s.point.r), site(s) {};
    };

    class CircleEvent : public Event {
    private:
        bool valid = true;
    public:
        // references to the beach line elements that cause the circle event
        rBeachLineElement first, second;
        // the point at which a Voronoi vertex is created when the event happens
        const Point center;
        CircleEvent(rBeachLineElement first, rBeachLineElement second, Point center, _float_t r)
                : Event(r + center.r), first(first), second(second), center(center) {};

        void invalidate() {
            valid = false;
        }
        [[nodiscard]] bool isValid() const {
            return valid;
        }
    };

    // comparison of events based on their radius
    template <class T>
    struct CmpEvent {
        static_assert(is_base_of<Event, T>::value, "template type T must inherit from Event");
        static bool greater_equal(const Event* e1, const Event* e2) {
            return e1->r > e2->r;
        }

        bool operator()(const shared_ptr<T>& e1, const shared_ptr<T> &e2) {
            return greater_equal(e1.get(), e2.get());
        }
    };

    // queue of events
    template <class T>
    class EventQueue : public priority_queue<shared_ptr<T>, vector<shared_ptr<T>>, CmpEvent<T>> {
    public:
        const T* getTop() {
            return (this->empty()) ? nullptr : this->top().get();
        }
    };

    // struct identifying a set of three sites based on their IDs. Used for caching the prediction of circle events
    struct SiteTriple {
        unsigned long long ID1, ID2, ID3;
        SiteTriple(rSite a, rSite b, rSite c) {
            vector<unsigned long long> ids = {a.ID, b.ID, c.ID};
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
    using SiteTripleMap = unordered_map<SiteTriple, _Point<_float_T>, HashSiteTriple, CmpSiteTriple>;


}