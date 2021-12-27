#include <iostream>

#include <lib.hpp>
#include <datastructures.hpp>
#include <calculations.hpp>

namespace hyperbolic {
    _float_t BeachLineElement::calculateAngularCoordinate(_float_t r_sweep) {
        return calculate_beach_line_intersection(&first, &second, r_sweep);
    }

    int BeachLineElement::getCount(const BeachLineElement * const e) {
        return (e) ? e->count : 0;
    }

    void BeachLineElement::setParent(pBeachLineElement e, BeachLineElement * p) {
        if (e) e->parent = p;
    }

    void BeachLineElement::update(pBeachLineElement e) {
        if (e) {
            e->count = 1 + getCount(e->leftChild) + getCount(e->rightChild);
            setParent(e->leftChild, e);
            setParent(e->rightChild, e);
        }
    }

    void BeachLineElement::split(pBeachLineElement t, pBeachLineElement& l, pBeachLineElement& r, int pos, int add=0) {
        // splits such that l holds elements from positions [0 .. pos-1], and r from [pos .. n]
        if (!t) {
            l = r = nullptr;
            return;
        }

        int current_position = getCount(t->leftChild) + add;
        if (pos > current_position) {
            split(t->rightChild, t->rightChild, r, pos, current_position + 1);
            l = t;
        } else if (pos > current_position){
            split(t->leftChild, l, t->leftChild, pos, add);
            r = t;
        } else {
            l = t->leftChild;
            t->leftChild = nullptr;
            r = t;
        }
        update(t);
    }

    void BeachLineElement::merge(pBeachLineElement& t, pBeachLineElement l, pBeachLineElement r) {
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

    void BeachLineElement::find(pBeachLineElement& result, pBeachLineElement t, int pos, int add=0) {
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

    void BeachLineElement::find(pBeachLineElement& result, int& pos, pBeachLineElement t, _float_t theta, _float_t r_sweep, int add=0) {
        // returns the first element clockwise of theta in t (i.e. with a angular coordinate smaller/equal theta)
        if (!t) {
            result = nullptr;
            return;
        }
        int current_position = getCount(t->leftChild) + add;
        _float_t current_key = t->calculateAngularCoordinate(r_sweep);
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

    int BeachLineElement::getPosition(pBeachLineElement e) {
        int pos = getCount(e->leftChild);
        while (e->parent) {
            e = e->parent;
            pos += getCount(e->leftChild) + 1;
        }
        return pos;
    }

    void BeachLineElement::del(pBeachLineElement t) {
        if (!t) return;
        del(t->leftChild);
        del(t->rightChild);
        delete t;
    }

    BeachLinePosition BeachLine::search(const SiteEvent & s) {
        pBeachLineElement first; int position_first;
        BeachLineElement::find(first, position_first, root, s.site.theta, s.r);

        auto size = BeachLineElement::getCount(root);
        if (!first) {
            BeachLineElement::find(first, root, 0);
            position_first = 0;
        }

        pBeachLineElement second;
        BeachLineElement::find(second, root, (position_first+1) % size);

        return BeachLinePosition(first, second, position_first);
    }

    void BeachLine::insert(int position, rBeachLineElement firstNew, rBeachLineElement secondNew) {
        // inserts the beach line elements firstNew and second new after position positionFirst and rearranges the order such that firstNew is the first element and secondNew is the last element
        auto size = BeachLineElement::getCount(root);

        pBeachLineElement left, right;
        BeachLineElement::split(root, left, right, (position + 1) % size);

        BeachLineElement::merge(right, &firstNew, right);
        BeachLineElement::merge(left, left, &secondNew);
        BeachLineElement::merge(root, right, left);
    }

    void BeachLine::replace(const CircleEvent & e, rBeachLineElement newElement) {
        // replaces the elements specified in e and replaces them with newElement
        // Note that we can safely assume that e does not refer to the last and first element as our arrangement can never have these two involved in a circle event

        int position = BeachLineElement::getPosition(&e.first);
        pBeachLineElement left, middle, right;
        BeachLineElement::split(root, left, middle, position);
        BeachLineElement::split(middle, middle, right, 2);

        delete middle;

        BeachLineElement::merge(left, left, &newElement);
        BeachLineElement::merge(root, left, right);
    }

    int BeachLine::size() {
        return BeachLineElement::getCount(root);
    }

#ifndef NDEBUG
    void BeachLine::print(pBeachLineElement t) {
        if (!t) return;
        print(t->leftChild);
        std::cout << "(" << t->first.ID  << ") ";
        print(t->rightChild);
        std::cout << std::endl;
    };
#endif
}