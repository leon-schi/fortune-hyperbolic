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

    typedef BeachLineElement * pBeachLineElement;
    typedef BeachLineElement & rBeachLineElement;

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
        CircleEvent(rBeachLineElement first, rBeachLineElement second, Point center, _float_t r) : Event(r + center.r), first(first), second(second), center(center) {};

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
    using SiteTripleMap = unordered_map<SiteTriple, Point, HashSiteTriple, CmpSiteTriple>;

    //TODO: ensure efficient access to neighbors

    /*
     * Elements that are within the beach line.
     * Represented as a tuple of sites (a,b) indicating that this element represents the point where the active segment changes from a to b in ccw direction
     * */
    class BeachLineElement {
        friend class BeachLine;
    public:
        // the two sites defining the element
        rSite first, second;
        // pointer to circle element of this element with the next / previous element in the beach line
        shared_ptr<CircleEvent> next = nullptr, previous = nullptr;
        // pointer to an edge on which the element moves
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

        // used to assign a vertex to the right vertex from the stored edge when this element is involved in a circle event
        void assignVertex(rPoint v) {assignableVertex = &v;}
    private:
        static void del(pBeachLineElement e);
        // pointer to a reference of one of the vertices int the stored edge
        pPoint& assignableVertex;

        // internal properties used for the data structure
        pBeachLineElement leftChild = nullptr, rightChild = nullptr;
        pBeachLineElement parent = nullptr;

        int count = 0;
        uint32_t priority;
    };

    /*
     * class implementing all operations on the beach line
     * currently implemented as a treap
     * */
    class BeachLine {
    private:
        // root of the treap
        pBeachLineElement root = nullptr;
        // angular coordinate of the last inserted vertex
        _float_t reference_angle = 0;

        // transforms an angle to the corresponding angle under the perspective of reference_angle
        [[nodiscard]] _float_t transformAngle(_float_t theta) const;

        // functions for updating properties of BeachLineElements
        static void setParent(pBeachLineElement e, pBeachLineElement p);
        [[nodiscard]] static int getCount(const BeachLineElement* e);
        static void update(pBeachLineElement e);

        // split and merge as fundamental operation of the treap
        static void split(pBeachLineElement t, pBeachLineElement& l, pBeachLineElement& r, int pos, int add);
        static void merge(pBeachLineElement& t, pBeachLineElement l, pBeachLineElement r);

        // returns the position of a given BeachLineELement
        static int getPosition(pBeachLineElement e);
        // copies a pointer to all remaining elements in the tree rooted at e into v
        static void getRemainingElements(pBeachLineElement e, vector<pBeachLineElement>& v);

        // wrapper for calculating the current angular coordinate of a BeachLineElement
        [[nodiscard]]  _float_t calculateAngularCoordinate(pBeachLineElement e, _float_t r_sweep, int pos) const;

        // search for the element at position pos in the tree rooted at t
        void find(pBeachLineElement& result, pBeachLineElement t, int pos, int add);
        // binary search returning the last BeachLineElement with an angular coordinate smaller than theta + its position at a sweep circle radius of r_sweep
        void find(pBeachLineElement& result, int& position, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add);

        // finds the leftmost/rightmost child of the tree rooted at e
        void getLeftmostChild(pBeachLineElement& result, pBeachLineElement e);
        void getRightmostChild(pBeachLineElement& result, pBeachLineElement e);
    public:
        BeachLine() = default;
        ~BeachLine() {
            delete root;
        }

        /*
         * returns the position of the first BeachLineElement clockwise of s.theta. first and second are set as
         * the first element clockwise and counterclockwise of s.theta, respectively
         * */
        int search(rPoint s, _float_t r_sweep, pBeachLineElement& first, pBeachLineElement& second);

        /*
         * inserts the elements firstNew and second new right after the position-th element.
         * It then re-arranges the beach line such that firstNew is the first element and
         * second new is the last element. firstNew.second must be equal to secondNew.first
         * */
        void insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew);

        /*
         * removes the elements specified in e and replaces them with newElement.
         * leftNeighbor and rightNeighbor are set to be the left (clockwise) and right (counterclockwise)
         * neighbor of the new element.
         * */
        void replace(const CircleEvent& e, rBeachLineElement newElement, pBeachLineElement& leftNeighbor, pBeachLineElement& rightNeighbor);
        [[nodiscard]] int size() const;

        /*
         * copies a pointer to all remaining elements into v
         * */
        void getRemainingElements(vector<pBeachLineElement>& elements);

#ifndef NDEBUG
        // methods for printing the beach line when debugging
        void print(_float_t r_sweep);
        void print(pBeachLineElement t, _float_t r_sweep, int add);
#endif
    };
}