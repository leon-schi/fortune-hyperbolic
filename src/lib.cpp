#include <iostream>
#include <queue>

#include <lib.hpp>
#include <fortune.hpp>
#include <datastructures.hpp>
#include <calculations.hpp>
#include <cmake.hpp>

using std::cout, std::endl, std::unique_ptr, std::make_unique, std::vector;
using namespace hyperbolic;

namespace hyperbolic {
    void print_version() {
        cout << VERSION_MAJOR << endl;
    }

    FortuneHyperbolicImplementation::FortuneHyperbolicImplementation(VoronoiDiagram &v, const vector<Point> &sites) : voronoiDiagram(v) {
        unsigned long long id = 0;
        for (const Point& p : sites) {
            this->sites.emplace_back(Site(p, id));
            id++;
        }
    }

    unique_ptr<FortuneHyperbolic> getInstance(VoronoiDiagram& diagram, const std::vector<Point>& sites) {
        return make_unique<FortuneHyperbolicImplementation>(diagram, sites);
    }

    void invalidateCircleEvent(CircleEvent* e) {
        if (e) e->invalidate();
    }

    void FortuneHyperbolicImplementation::addCircleEvent(rBeachLineElement first, rBeachLineElement second) {
        invalidateCircleEvent(first.next);
        invalidateCircleEvent(second.previous);
        Point centerCircleEvent;
        if (predict_circle_event(centerCircleEvent, first, second)) {
            _float_t radius = distance(centerCircleEvent, first.first.point);
            unique_ptr<CircleEvent> ce = make_unique<CircleEvent>(first, second, centerCircleEvent, radius);
            first.next = ce.get();
            second.previous = ce.get();
            circleEventQueue.push(std::move(ce));
        }
    }

    void FortuneHyperbolicImplementation::handleSiteEvent(const SiteEvent& e) {
#ifndef NDEBUG
        cout << "Handling Site Event with radius " << e.r << endl;
#endif

        BeachLinePosition p = beachLine.search(e);
        rSite hitSite = p.first->second;
        Edge& edge = getNewEdge(e.site, hitSite);
        auto first = new BeachLineElement(e.site, hitSite, edge, edge.firstVertex);
        auto last = new BeachLineElement(hitSite, e.site, edge, edge.secondVertex);

        addCircleEvent(*p.first, *last);
        addCircleEvent(*first, *p.second);

        beachLine.insert(p.positionFirst, *first, *last);
    }

    void FortuneHyperbolicImplementation::handleCircleEvent(const CircleEvent& e) {
        if (!e.isValid()) return;

#ifndef NDEBUG
        cout << "Handling Circle Event with radius " << e.r << endl;
#endif

        rPoint v = getNewVertex(e.center);
        e.first.assignVertex(v);
        e.second.assignVertex(v);

        // replace first and second with new beach line element
        rSite a = e.first.first, b = e.second.second;
        Edge& edge = getNewEdge(a, b);
        edge.firstVertex = &v;
        auto newElement = new BeachLineElement(a, b, edge, edge.secondVertex);

        pBeachLineElement leftNeighbor, rightNeighbor;
        beachLine.replace(e, *newElement, leftNeighbor, rightNeighbor);

        // predict new circle events
        addCircleEvent(*leftNeighbor, *newElement);
        addCircleEvent(*rightNeighbor, *newElement);
    }

    bool FortuneHyperbolicImplementation::eventsRemaining() {
        return !(siteEventQueue.empty() && circleEventQueue.empty());
    }

    void FortuneHyperbolicImplementation::initializeBeachLine() {
        rSite a = siteEventQueue.top()->site;
        siteEventQueue.pop();
        rSite b = siteEventQueue.top()->site;
        siteEventQueue.pop();

        Edge& edge = getNewEdge(a, b);
        auto first = new BeachLineElement(b, a, edge, edge.firstVertex);
        auto last = new BeachLineElement(a, b, edge, edge.secondVertex);
        beachLine.insert(0, *first, *last);
    }

    void FortuneHyperbolicImplementation::calculate() {
        if (isCalculated) return;
        isCalculated = true;

        if (sites.size() <= 1) return;

        for (rSite s : sites) {
            siteEventQueue.push(make_unique<SiteEvent>(s));
        }

        initializeBeachLine();

        while (eventsRemaining()) {
            const SiteEvent* se = siteEventQueue.getTop();
            const CircleEvent* ce = circleEventQueue.getTop();

            if (!ce || (se && se->r <= ce->r)) {
                handleSiteEvent(*se);
                siteEventQueue.pop();
            } else {
                handleCircleEvent(*ce);
                circleEventQueue.pop();
            }

#ifndef NDEBUG
            beachLine.print();
#endif
        }
    }

    Edge& FortuneHyperbolicImplementation::getNewEdge(rSite a, rSite b) {
        voronoiDiagram.edges.emplace_back(Edge(a.point, b.point));
        return voronoiDiagram.edges.back();
    }

    Point& FortuneHyperbolicImplementation::getNewVertex(rPoint p) {
        voronoiDiagram.vertices.emplace_back(p);
        return voronoiDiagram.vertices.back();
    }
}