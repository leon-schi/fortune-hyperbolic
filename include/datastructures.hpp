#pragma once

#include <vector>
#include <queue>
#include <random>
#include <type_traits>

#include <lib.hpp>
#include <calculations.hpp>

using std::is_base_of;

namespace hyperbolic {
    class BeachLineElement;

    typedef BeachLineElement * pBeachLineElement;
    typedef BeachLineElement & rBeachLineElement;

    class Event {
    public:
        const _float_t r {0};
        explicit Event(_float_t r) : r(r) {};
        virtual ~Event() = default;
        [[nodiscard]] _float_t getRadius() const {
            return r;
        }
    };

    class SiteEvent : public Event {
    public:
        rSite site;
        explicit SiteEvent(rSite s) : Event(s.r), site(s) {};
    };

    class CircleEvent : public Event {
    private:
        bool valid = true;
    public:
        rBeachLineElement first, second;
        CircleEvent(rBeachLineElement first, rBeachLineElement second, _float_t r) : Event(r), first(first), second(second) {};
        void invalidate() {
            valid = false;
        }
        [[nodiscard]] bool isValid() const {
            return valid;
        }
    };

    template <class T>
    struct CmpEvent {
        static_assert(is_base_of<Event, T>::value, "template type T must inherit from Event");
        bool compare(const Event* e1, const Event* e2) {
            return e1->r > e2->r;
        }

        bool operator()(const std::unique_ptr<T> &e1, const std::unique_ptr<T> &e2) {
            return compare(e1.get(), e2.get());
        }
    };

    template <class T>
    class EventQueue : public std::priority_queue<std::unique_ptr<T>, std::vector<std::unique_ptr<T>>, CmpEvent<T>> {
    public:
        const T* getTop() {
            return (this->empty()) ? nullptr : this->top().get();
        }
    };

    struct BeachLinePosition {
        pBeachLineElement first, second;
        int positionFirst;
        BeachLinePosition(pBeachLineElement first, pBeachLineElement second, int positionFirst) :
            first(first), second(second), positionFirst(positionFirst) {}
    };

    class BeachLineElement {
        friend class BeachLine;
    private:
        rSite first, second;
        const Edge & edge;
        pVertex& assignableVertex;
        CircleEvent *left = nullptr, *right = nullptr;

        pBeachLineElement leftChild = nullptr, rightChild = nullptr;
        pBeachLineElement parent = nullptr;

        int count = 0;
        uint32_t priority;

        [[nodiscard]] _float_t calculateAngularCoordinate(_float_t r_sweep);

        static void setParent(pBeachLineElement e, pBeachLineElement p);
        [[nodiscard]] static int getCount(const BeachLineElement* e);
        static void update(pBeachLineElement e);
        static void split(pBeachLineElement t, pBeachLineElement& l, pBeachLineElement& r, int pos, int add);
        static void merge(pBeachLineElement& t, pBeachLineElement l, pBeachLineElement r);
        static void find(pBeachLineElement& result, pBeachLineElement t, int pos, int add);
        static void find(pBeachLineElement& result, int& position, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add);
        static int getPosition(pBeachLineElement e);

        static void del(pBeachLineElement t);
    public:
        BeachLineElement(rSite first, rSite second, Edge& edge, pVertex& v) : first(first), second(second), edge(edge), assignableVertex(v) {
            std::random_device d;
            std::mt19937 rng(d());
            std::uniform_int_distribution<uint32_t> uint_dist;
            priority = uint_dist(rng);
        };
        ~BeachLineElement() {
            del(this);
        }
    };

    class BeachLine {
    private:
        pBeachLineElement root = nullptr;
    public:
        BeachLine() = default;
        BeachLinePosition search(const SiteEvent & s);
        void insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew);
        void replace(const CircleEvent & e, rBeachLineElement newElement);
        int size();

#ifndef NDEBUG
        void print(pBeachLineElement t);
#endif
    };
}