#include <iostream>

#include <lib.hpp>
#include <datastructures.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    _float_t BeachLine::calculateAngularCoordinate(pBeachLineElement e, _float_t r_sweep, int pos) const {
        if (pos == 0 && r_sweep == e->first.point.r)
            return 0;
        else if (pos == size() - 1 && r_sweep == e->second.point.r)
            return 2*M_PI;

        _float_t theta = calculate_beach_line_intersection(&e->first, &e->second, r_sweep);
        return clip(theta - reference_angle);
    }

    int BeachLine::getCount(const BeachLineElement * const e) {
        return (e) ? e->count : 0;
    }

    void BeachLine::setParent(pBeachLineElement e, BeachLineElement * p) {
        if (e) e->parent = p;
    }

    void BeachLine::update(pBeachLineElement e) {
        if (e) {
            e->count = 1 + getCount(e->leftChild) + getCount(e->rightChild);
            setParent(e->leftChild, e);
            setParent(e->rightChild, e);
            e->parent = nullptr;
        }
    }

    void BeachLine::split(pBeachLineElement t, pBeachLineElement& l, pBeachLineElement& r, int pos, int add=0) {
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
    }

    void BeachLine::merge(pBeachLineElement& t, pBeachLineElement l, pBeachLineElement r) {
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

    void BeachLine::find(pBeachLineElement& result, pBeachLineElement t, int pos, int add=0) {
        if (!t) {
            result = nullptr;
            return;
        }

        int current_position = getCount(t->leftChild) + add;
        if (pos == current_position)
            result = t;
        else if (pos > current_position)
            find(result, t->rightChild, pos, current_position+1);
        else
            find(result, t->leftChild, pos, add);
    }

    void BeachLine::find(pBeachLineElement& result, int& pos, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add=0) {
        // returns the first element clockwise of theta in t (i.e. with a angular coordinate smaller/equal theta)
        if (!t) {
            result = nullptr;
            return;
        }
        int current_position = getCount(t->leftChild) + add;
        _float_t current_key = calculateAngularCoordinate(t, r_sweep, current_position);
        if (theta >= current_key) {
            find(result, pos, t->rightChild, theta, r_sweep, current_position+1);
            if (result == nullptr) {
                result = t;
                pos = current_position;
            }
            return;
        } else {
            find(result, pos, t->leftChild, theta, r_sweep, add);
        }
    }

    int BeachLine::getPosition(pBeachLineElement e) {
        int pos = getCount(e->leftChild);
        while (e->parent) {
            if (e->parent->leftChild != e)
                pos += getCount(e->parent->leftChild) + 1;
            e = e->parent;
        }
        return pos;
    }

    void BeachLineElement::del(pBeachLineElement t) {
        if (!t) return;
        del(t->leftChild);
        t->leftChild = nullptr;
        del(t->rightChild);
        t->rightChild = nullptr;
        delete t;
    }

    void BeachLine::getLeftmostChild(pBeachLineElement& result, pBeachLineElement e) {
        if (!e) {
            result = nullptr;
            return;
        }
        find(result, e, 0);
    };

    void BeachLine::getRightmostChild(pBeachLineElement &result, pBeachLineElement e) {
        if (!e) {
            result = nullptr;
            return;
        }
        find(result, e, getCount(e)-1);
    };

    BeachLinePosition BeachLine::search(const SiteEvent & s) {
        pBeachLineElement first; int position_first;
        _float_t theta = clip(s.site.point.theta - reference_angle);
        BeachLine::find(first, position_first, root, theta, s.r);

        auto size = BeachLine::getCount(root);
        if (!first) {
            BeachLine::find(first, root, 0);
            position_first = 0;
        }

        pBeachLineElement second;
        BeachLine::find(second, root, (position_first+1) % size);

        return BeachLinePosition(first, second, position_first);
    }

    void BeachLine::insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew) {
        // inserts the neighboring beach line elements firstNew and second new after position positionFirst and rearranges the order such that firstNew is the first element and secondNew is the last element
        auto size = BeachLine::getCount(root);
        position = (size > 0) ? (position + 1) % size : 0;

        pBeachLineElement left, right;
        BeachLine::split(root, left, right, position);

        BeachLine::merge(right, &firstNew, right);
        BeachLine::merge(left, left, &secondNew);
        BeachLine::merge(root, right, left);

        reference_angle = firstNew.first.point.theta;
    }

    void BeachLine::replace(const CircleEvent & e, rBeachLineElement newElement, pBeachLineElement& leftNeighbor, pBeachLineElement& rightNeighbor) {
        // replaces the elements specified in e and replaces them with newElement
        // Note that we can safely assume that e does not refer to the last and first element as our arrangement can never have these two involved in a circle event

        int position = BeachLine::getPosition(&e.first);
        pBeachLineElement left, middle, right;
        BeachLine::split(root, left, middle, position);
        BeachLine::split(middle, middle, right, 2);

        delete middle;

        // get left and right neighbors
        BeachLine::getRightmostChild(leftNeighbor, left);
        BeachLine::getLeftmostChild(rightNeighbor, right);
        if (leftNeighbor == nullptr) BeachLine::getRightmostChild(leftNeighbor, right);
        if (rightNeighbor == nullptr) BeachLine::getLeftmostChild(rightNeighbor, left);

        BeachLine::merge(left, left, &newElement);
        BeachLine::merge(root, left, right);
    }

    int BeachLine::size() const {
        return BeachLine::getCount(root);
    }

    void BeachLine::getRemainingElements(pBeachLineElement e, vector<pBeachLineElement>& v) {
        if (!e) return;
        getRemainingElements(e->leftChild, v);
        v.push_back(e);
        getRemainingElements(e->rightChild, v);
    }

    void BeachLine::getRemainingElements(vector<pBeachLineElement>& v) {
        getRemainingElements(root, v);
    }

#ifndef NDEBUG
    void BeachLine::print(_float_t r_sweep) {
        print(root, r_sweep, 0);
        std::cout << std::endl;
    }

    void BeachLine::print(pBeachLineElement t, _float_t r_sweep, int add=0) {
        if (!t) return;
        print(t->leftChild, r_sweep, add);
        int cur_pos = getCount(t->leftChild) + add;
        std::cout << "((" << t->first.ID << ", " << t->second.ID << "), " << cur_pos << ", " << calculateAngularCoordinate(t, r_sweep, cur_pos) << ") ";
        print(t->rightChild, r_sweep, cur_pos + 1);
    };
#endif
}