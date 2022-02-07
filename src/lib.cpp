#include <iostream>
#include <queue>
#include <fstream>

#include <lib.hpp>
#include <fortune.hpp>
#include <datastructures.hpp>
#include <calculations.hpp>
#include <cmake.hpp>

using std::cout, std::endl, std::make_shared, std::make_unique, std::vector;
using namespace hyperbolic;

namespace hyperbolic {
    FortuneHyperbolicImplementation::FortuneHyperbolicImplementation(VoronoiDiagram &v, const vector<Point> &sites) : r_sweep(0), voronoiDiagram(v) {
        unsigned long long id = 0;
        for (const Point& p : sites) {
            this->sites.emplace_back(Site(p, id));
            id++;
        }
    }

    unique_ptr<FortuneHyperbolic> getNewInstance(VoronoiDiagram& diagram, const std::vector<Point>& sites) {
        return make_unique<FortuneHyperbolicImplementation>(diagram, sites);
    }

    void invalidateCircleEvent(CircleEvent* e) {
        if (e) e->invalidate();
    }

    void FortuneHyperbolicImplementation::addCircleEvent(rBeachLineElement first, rBeachLineElement second) {
        invalidateCircleEvent(first.next.get());
        invalidateCircleEvent(second.previous.get());
        Point centerCircleEvent;
        if (predictCircleEvent(centerCircleEvent, first, second)) {
            _float_t radius = distance(centerCircleEvent, first.first.point);
            shared_ptr<CircleEvent> ce = make_shared<CircleEvent>(first, second, centerCircleEvent, radius);
            if (ce->r >= r_sweep) {
                first.next = ce;
                second.previous = ce;
                circleEventQueue.push(ce);
            }
        }
    }

    void FortuneHyperbolicImplementation::handleSiteEvent(const SiteEvent& e) {
#ifndef NDEBUG
        cout << "Handling Site Event with radius " << e.r << endl;
#endif
        r_sweep = e.r;

        pBeachLineElement first, second;
        int positionFirst = beachLine.search(e.site.point, e.r, first, second);
        rSite hitSite = first->second;

        Edge* edge = getNewEdge(e.site, hitSite, EdgeType::BIDIRECTIONAL);
        auto firstNew = new BeachLineElement(e.site, hitSite, edge, edge->firstVertex);
        auto secondNew = new BeachLineElement(hitSite, e.site, edge, edge->secondVertex);

        addCircleEvent(*first, *secondNew);
        addCircleEvent(*firstNew, *second);

        beachLine.insert(positionFirst, *firstNew, *secondNew);
    }

    void FortuneHyperbolicImplementation::handleCircleEvent(const CircleEvent& e) {
#ifndef NDEBUG
        if (!e.isValid())
            cout << "Skipping invalid circle event" << endl;
        else
            cout << "Handling Circle Event with radius " << e.r << " concerning beach line elements (" << e.first.first.ID << ", " << e.first.second.ID << ") and (" << e.second.first.ID << ", " << e.second.second.ID << ")" << endl;
#endif

        if (!e.isValid()) return;
        r_sweep = e.r;

        pPoint v = getNewVertex(e.center);

        e.first.assignVertex(*v);
        e.second.assignVertex(*v);

        // replace first and second with new beach line element
        rSite a = e.first.first, b = e.second.second;
        Edge* edge = getNewEdge(a, b, (a.point.r >= b.point.r) ? EdgeType::CCW : EdgeType::CW);
        edge->firstVertex = v;
        auto newElement = new BeachLineElement(a, b, edge, edge->secondVertex);

        pBeachLineElement leftNeighbor, rightNeighbor;
        beachLine.replace(e, *newElement, leftNeighbor, rightNeighbor);

        // predict new circle events
        addCircleEvent(*leftNeighbor, *newElement);
        addCircleEvent(*newElement, *rightNeighbor);
    }

    bool FortuneHyperbolicImplementation::eventsRemaining() {
        return !(siteEventQueue.empty() && circleEventQueue.empty());
    }

    void FortuneHyperbolicImplementation::initializeBeachLine() {
        // TODO: handle speicial case that the tree closest sites are equidistant from the origin

        rSite a = siteEventQueue.top()->site;
        siteEventQueue.pop();
        rSite b = siteEventQueue.top()->site;
        siteEventQueue.pop();
        r_sweep = b.point.r;

        Edge* edge = getNewEdge(b, a, EdgeType::BIDIRECTIONAL);
        auto first = new BeachLineElement(b, a, edge, edge->firstVertex);
        auto last = new BeachLineElement(a, b, edge, edge->secondVertex);
        beachLine.insert(0, *first, *last);
    }

    void FortuneHyperbolicImplementation::calculate() {
        if (isCalculated) return;
        isCalculated = true;

        if (sites.size() <= 1) return;

        for (rSite s : sites) {
            siteEventQueue.push(make_shared<SiteEvent>(s));
        }

        initializeBeachLine();
#ifndef NDEBUG
        beachLine.print(r_sweep);
#endif

        while (eventsRemaining()) {
            const SiteEvent* se = siteEventQueue.getTop();
            const CircleEvent* ce = circleEventQueue.getTop();

            if (!ce || (se && CmpEvent<Event>::greater_equal(ce, se))) {
                handleSiteEvent(*se);
                siteEventQueue.pop();
            } else {
                handleCircleEvent(*ce);
                circleEventQueue.pop();
            }

#ifndef NDEBUG
            beachLine.print(r_sweep);
#endif
        }

        sanitizeEdges();
    }

    void FortuneHyperbolicImplementation::sanitizeEdges() {
        vector<pBeachLineElement> remaining;
        beachLine.getRemainingElements(remaining);
        for (pBeachLineElement e : remaining) {
            _float_t theta = calculate_beach_line_intersection(&e->first, &e->second, r_sweep + 1);
            Bisector b(&e->first.point, &e->second.point);
            e->edge->last_known_position = Point(b(theta), theta);
        }
    }

    Edge* FortuneHyperbolicImplementation::getNewEdge(rSite a, rSite b, EdgeType edgeType) {
        voronoiDiagram.edges.emplace_back(make_unique<Edge>(a, b, edgeType));
        return voronoiDiagram.edges.back().get();
    }

    Point* FortuneHyperbolicImplementation::getNewVertex(rPoint p) {
        voronoiDiagram.vertices.push_back(make_unique<Point>(p));
        return voronoiDiagram.vertices.back().get();
    }
}