#pragma once
#include <iostream>
#include <datastructures.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    //TODO: ensure efficient access to neighbors

    /*
     * Elements that are within the beach line.
     * Represented as a tuple of sites (a,b) indicating that this element represents the point where the active segment changes from a to b in ccw direction
     * */
    class BeachLineElement {
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
        void assignVertex(rPoint v) const {
            assignableVertex = &v;
        }
        static void del(pBeachLineElement t) {
            if (!t) return;
            del(t->leftChild);
            t->leftChild = nullptr;
            del(t->rightChild);
            t->rightChild = nullptr;
            delete t;
        };
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
    template<class K>
    class BeachLine {
    private:
        K& kernel;
        // root of the treap
        pBeachLineElement root = nullptr;
        // angular coordinate of the last inserted vertex
        _float_t reference_angle = 0;

        // transforms an angle to the corresponding angle under the perspective of reference_angle
        [[nodiscard]] _float_t transformAngle(_float_t theta) const {
            return clip<_float_t>(theta - reference_angle);
        };

        // functions for updating properties of BeachLineElements
        static void setParent(pBeachLineElement e, pBeachLineElement p) {
            if (e) e->parent = p;
        };

        [[nodiscard]] static int getCount(pBeachLineElement e) {
            return (e) ? e->count : 0;
        };

        static void update(pBeachLineElement e) {
            if (e) {
                e->count = 1 + getCount(e->leftChild) + getCount(e->rightChild);
                setParent(e->leftChild, e);
                setParent(e->rightChild, e);
                e->parent = nullptr;
            }
        };

        // split and merge as fundamental operation of the treap
        static void split(pBeachLineElement t, pBeachLineElement &l, pBeachLineElement &r, int pos, int add=0) {
            // splits such that l holds elements from positions [0 .. pos-1], and r from [pos .. n]
            if (!t) {
                l = r = nullptr;
                return;
            }

            int current_position = getCount(t->leftChild) + add;
            if (pos > current_position) {
                split(t->rightChild, t->rightChild, r, pos, current_position + 1);
                l = t;
            } else if (pos < current_position) {
                split(t->leftChild, l, t->leftChild, pos, add);
                r = t;
            } else {
                l = t->leftChild;
                t->leftChild = nullptr;
                r = t;
            }
            update(t);
        };

        static void merge(pBeachLineElement &t, pBeachLineElement l, pBeachLineElement r) {
            if (!l || !r)
                t = l ? l : r;
            else if (l->priority > r->priority) {
                merge(l->rightChild, l->rightChild, r);
                t = l;
            } else {
                merge(r->leftChild, l, r->leftChild);
                t = r;
            }

            update(t);
        };

        // returns the position of a given BeachLineElement
        static int getPosition(pBeachLineElement e) {
            int pos = getCount(e->leftChild);
            while (e->parent) {
                if (e->parent->leftChild != e)
                    pos += getCount(e->parent->leftChild) + 1;
                e = e->parent;
            }
            return pos;
        };

        // copies a pointer to all remaining elements in the tree rooted at e into v
        static void getRemainingElements(pBeachLineElement e, vector<pBeachLineElement> &v) {
            if (!e) return;
            getRemainingElements(e->leftChild, v);
            v.push_back(e);
            getRemainingElements(e->rightChild, v);
        };

        // search for the element at position pos in the tree rooted at t
        void find(pBeachLineElement &result, pBeachLineElement t, int pos, int add=0) {
            if (!t) {
                result = nullptr;
                return;
            }

            int current_position = getCount(t->leftChild) + add;
            if (pos == current_position)
                result = t;
            else if (pos > current_position)
                find(result, t->rightChild, pos, current_position + 1);
            else
                find(result, t->leftChild, pos, add);
        };

        // binary search returning the last BeachLineElement with an angular coordinate smaller than theta + its position at a sweep circle radius of r_sweep
        void find(pBeachLineElement &result, int &pos, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add=0) {
            // returns the first element clockwise of theta in t (i.e. with a angular coordinate smaller/equal theta)
            if (!t) {
                result = nullptr;
                return;
            }
            int current_position = getCount(t->leftChild) + add;

            if (kernel.before(theta, reference_angle, t->first.point, t->second.point, r_sweep)) {
                find(result, pos, t->leftChild, theta, r_sweep, add);
            } else {
                find(result, pos, t->rightChild, theta, r_sweep, current_position + 1);
                if (result == nullptr) {
                    result = t;
                    pos = current_position;
                }
                return;
            }
        };

        // finds the leftmost/rightmost child of the tree rooted at e
        void getLeftmostChild(pBeachLineElement &result, pBeachLineElement e) {
            if (!e) {
                result = nullptr;
                return;
            }
            find(result, e, 0);
        };

        void getRightmostChild(pBeachLineElement &result, pBeachLineElement e) {
            if (!e) {
                result = nullptr;
                return;
            }
            find(result, e, getCount(e) - 1);
        };

        _float_t calculateAngularCoordinate(pBeachLineElement e, _float_t r_sweep, int pos) const {
            // if e is the first or last element and if the sweep circle is very close to one of its sites, we have to be careful
            rPoint outer = (e->first.point.r >= e->second.point.r) ? e->first.point : e->second.point;
            if (r_sweep - outer.r <= 1e-20) {
                if (pos == 0) return 0;
                else if (pos == size() -1) return 2*M_PI;
                else return outer.theta;
            }

            _float_t theta = calculate_beach_line_intersection(&e->first, &e->second, r_sweep);
            return transformAngle(theta);
        }
    public:
        BeachLine(K& k) : kernel(k) {};

        ~BeachLine() {
            delete root;
        }

        /*
         * returns the position of the first BeachLineElement clockwise of s.theta. first and second are set as
         * the first element clockwise and counterclockwise of s.theta, respectively
         * */
        int search(rPoint s, _float_t r_sweep, pBeachLineElement &first, pBeachLineElement &second) {
            first = nullptr;
            int position_first;
            find(first, position_first, root, transformAngle(s.theta), r_sweep);

            auto size = BeachLine::getCount(root);
            if (!first) {
                position_first = getCount(root) - 1;
                BeachLine::find(first, root, position_first);
            }

            second = nullptr;
            BeachLine::find(second, root, (position_first + 1) % size);

            return position_first;
        };

        /*
         * inserts the elements firstNew and second new right after the position-th element.
         * It then re-arranges the beach line such that firstNew is the first element and
         * second new is the last element. firstNew.second must be equal to secondNew.first
         * */
        void insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew) {
            // inserts the neighboring beach line elements firstNew and second new after position positionFirst and rearranges the order such that firstNew is the first element and secondNew is the last element
            auto size = BeachLine::getCount(root);
            position = (size > 0) ? (position + 1) % size : 0;

            pBeachLineElement left, right;
            split(root, left, right, position);

            BeachLine::merge(right, &firstNew, right);
            BeachLine::merge(left, left, &secondNew);
            BeachLine::merge(root, right, left);

            reference_angle = firstNew.first.point.theta;
        };

        /*
         * removes the elements specified in e and replaces them with newElement.
         * leftNeighbor and rightNeighbor are set to be the left (clockwise) and right (counterclockwise)
         * neighbor of the new element.
         * */
        void replace(const CircleEvent &e, rBeachLineElement newElement, pBeachLineElement &leftNeighbor,
                     pBeachLineElement &rightNeighbor) {
            // replaces the elements specified in e and replaces them with newElement
            // Note that we can safely assume that e does not refer to the last and first element as our arrangement can never have these two involved in a circle event

            int position = BeachLine::getPosition(&e.first);
            pBeachLineElement left, middle, right;
            split(root, left, middle, position);
            split(middle, middle, right, 2);

            delete middle;

            // get left and right neighbors
            BeachLine::getRightmostChild(leftNeighbor, left);
            BeachLine::getLeftmostChild(rightNeighbor, right);
            if (leftNeighbor == nullptr) BeachLine::getRightmostChild(leftNeighbor, right);
            if (rightNeighbor == nullptr) BeachLine::getLeftmostChild(rightNeighbor, left);

            BeachLine::merge(left, left, &newElement);
            BeachLine::merge(root, left, right);
        };

        [[nodiscard]] int size() const {
            return BeachLine::getCount(root);
        };

#ifndef NDEBUG

        // methods for printing the beach line when debugging
        void print(_float_t r_sweep) {
            print(root, r_sweep, 0);
            std::cout << std::endl;
        };

        void print(pBeachLineElement t, _float_t r_sweep, int add) {
            if (!t) return;
            print(t->leftChild, r_sweep, add);
            int cur_pos = getCount(t->leftChild) + add;
            std::cout << "((" << t->first.ID << ", " << t->second.ID << "), " << cur_pos << ", "
                      << calculateAngularCoordinate(t, r_sweep, cur_pos) << ") ";
            print(t->rightChild, r_sweep, cur_pos + 1);
        };
#endif
    };
}