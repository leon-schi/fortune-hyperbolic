#pragma once

#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <type_traits>

#include <lib.hpp>

using std::is_base_of, std::unique_ptr, std::shared_ptr, std::vector, std::priority_queue, std::size_t, std::string, std::to_string, std::hash, std::unordered_map;

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

        bool operator()(const shared_ptr<T>& e1, const shared_ptr<T> &e2) {
            return compare(e1.get(), e2.get());
        }
    };

    template <class T>
    class EventQueue : public priority_queue<shared_ptr<T>, vector<shared_ptr<T>>, CmpEvent<T>> {
    public:
        const T* getTop() {
            return (this->empty()) ? nullptr : this->top().get();
        }
    };

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
    using SiteTripleMap = unordered_map<SiteTriple, Point, HashSiteTriple, CmpSiteTriple>;

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
        shared_ptr<CircleEvent> next = nullptr, previous = nullptr;
        Edge* const edge;
        BeachLineElement(rSite first, rSite second, Edge* edge, pPoint& v) : first(first), second(second), edge(edge), assignableVertex(v) {
            std::random_device d;
            std::mt19937 rng(d());
            //rng.seed(14232645);
            std::uniform_int_distribution<uint32_t> uint_dist;
            priority = uint_dist(rng);
        };
        ~BeachLineElement() {
            del(leftChild);
            del(rightChild);
        }
        void assignVertex(rPoint v) {assignableVertex = &v;}
    private:
        static void del(pBeachLineElement e);
        pPoint& assignableVertex;

        pBeachLineElement leftChild = nullptr, rightChild = nullptr;
        pBeachLineElement parent = nullptr;

        int count = 0;
        uint32_t priority;
    };

    class BeachLine {
    private:
        pBeachLineElement root = nullptr;
        _float_t reference_angle = 0; // angular coordinate of the reference point used for searching on the beach line

        static void setParent(pBeachLineElement e, pBeachLineElement p);
        [[nodiscard]] static int getCount(const BeachLineElement* e);
        static void update(pBeachLineElement e);

        static void split(pBeachLineElement t, pBeachLineElement& l, pBeachLineElement& r, int pos, int add);
        static void merge(pBeachLineElement& t, pBeachLineElement l, pBeachLineElement r);

        static int getPosition(pBeachLineElement e);

        static void getRemainingElements(pBeachLineElement e, vector<pBeachLineElement>& v);

        [[nodiscard]]  _float_t calculateAngularCoordinate(pBeachLineElement e, _float_t r_sweep, int pos) const;
        // search by known position in beach line
        void find(pBeachLineElement& result, pBeachLineElement t, int pos, int add);
        // binary search for handling site events
        void find(pBeachLineElement& result, int& position, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add);

        void getLeftmostChild(pBeachLineElement& result, pBeachLineElement e);
        void getRightmostChild(pBeachLineElement& result, pBeachLineElement e);
    public:
        BeachLine() = default;
        ~BeachLine() {
            delete root;
        }
        BeachLinePosition search(const SiteEvent & s);
        void insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew);
        void replace(const CircleEvent & e, rBeachLineElement newElement, pBeachLineElement& leftNeighbor, pBeachLineElement& rightNeighbor);
        [[nodiscard]] int size() const;

        void getRemainingElements(vector<pBeachLineElement>& elements);

#ifndef NDEBUG
        void print(_float_t r_sweep);
        void print(pBeachLineElement t, _float_t r_sweep, int add);
#endif
    };
}