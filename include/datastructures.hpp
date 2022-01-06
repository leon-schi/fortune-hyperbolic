#pragma once

#include <vector>
#include <queue>
#include <random>
#include <type_traits>

#include <lib.hpp>

using std::is_base_of, std::unique_ptr, std::vector, std::priority_queue;

namespace hyperbolic {
    class BeachLineElement;

    typedef BeachLineElement * pBeachLineElement;
    typedef BeachLineElement & rBeachLineElement;

    struct Site {
        rPoint point;
        const unsigned long long ID;
        Site(rPoint point, unsigned long long id) : point(point), ID(id) {};
    };
    typedef const Site & rSite;
    typedef const Site * pSite;

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
        rBeachLineElement first, second;
        const Point center;
        CircleEvent(rBeachLineElement first, rBeachLineElement second, Point center, _float_t r) : Event(r + center.r), first(first), second(second), center(center) {};
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

        bool operator()(const unique_ptr<T>& e1, const unique_ptr<T> &e2) {
            return compare(e1.get(), e2.get());
        }
    };

    template <class T>
    class EventQueue : public priority_queue<unique_ptr<T>, vector<unique_ptr<T>>, CmpEvent<T>> {
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

    //TODO: ensure efficient access to neighbors
    class BeachLineElement {
        friend class BeachLine;
    public:
        rSite first, second;
        CircleEvent *next = nullptr, *previous = nullptr;
        BeachLineElement(rSite first, rSite second, Edge& edge, pPoint& v) : first(first), second(second), edge(edge), assignableVertex(v) {
            std::random_device d;
            std::mt19937 rng(d());
            std::uniform_int_distribution<uint32_t> uint_dist;
            priority = uint_dist(rng);
        };
        ~BeachLineElement() {
            del(leftChild);
            del(rightChild);
        }
        void assignVertex(rPoint v) {assignableVertex = &v;}
    private:
        const Edge & edge;
        pPoint& assignableVertex;

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
        // search by known position in beach line
        static void find(pBeachLineElement& result, pBeachLineElement t, int pos, int add);
        // binary search for handling site events
        static void find(pBeachLineElement& result, int& position, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add);
        static int getPosition(pBeachLineElement e);
        static void del(pBeachLineElement);

        static void getLeftmostChild(pBeachLineElement& result, pBeachLineElement e);
        static void getRightmostChild(pBeachLineElement& result, pBeachLineElement e);
    };

    class BeachLine {
    private:
        pBeachLineElement root = nullptr;
    public:
        BeachLine() = default;
        ~BeachLine() {
            delete root;
        }
        BeachLinePosition search(const SiteEvent & s);
        void insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew);
        void replace(const CircleEvent & e, rBeachLineElement newElement, pBeachLineElement& leftNeighbor, pBeachLineElement& rightNeighbor);
        int size();

#ifndef NDEBUG
        void print();
        void print(pBeachLineElement t);
#endif
    };
}