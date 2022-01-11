#include <iostream>
#include <queue>

#include <lib.hpp>
#include <fortune.hpp>
#include <datastructures.hpp>
#include <calculations.hpp>
#include <cmake.hpp>

using std::cout, std::endl, std::make_shared, std::make_unique, std::vector;
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
        if (predict_circle_event(centerCircleEvent, first, second)) {
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

        BeachLinePosition p = beachLine.search(e);
        rSite hitSite = p.first->second;

        Bisector b(&e.site.point, &hitSite.point);

        Edge* edge = getNewEdge(e.site, hitSite);
        auto first = new BeachLineElement(e.site, hitSite, edge, edge->firstVertex);
        auto last = new BeachLineElement(hitSite, e.site, edge, edge->secondVertex);

        addCircleEvent(*p.first, *last);
        addCircleEvent(*first, *p.second);

        beachLine.insert(p.positionFirst, *first, *last);
    }

    void FortuneHyperbolicImplementation::handleCircleEvent(const CircleEvent& e) {
        if (!e.isValid()) return;
        r_sweep = e.r;

#ifndef NDEBUG
        cout << "Handling Circle Event with radius " << e.r << " concerning beach line elements (" << e.first.first.ID << ", " << e.first.second.ID << ") and (" << e.second.first.ID << ", " << e.second.second.ID << ")" << endl;
#endif

        pPoint v = getNewVertex(e.center);

        e.first.assignVertex(*v);
        e.second.assignVertex(*v);

        // replace first and second with new beach line element
        rSite a = e.first.first, b = e.second.second;
        Edge* edge = getNewEdge(a, b);
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
        
        Edge* edge = getNewEdge(a, b);
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

            if (!ce || (se && se->r <= ce->r)) {
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

    Edge* FortuneHyperbolicImplementation::getNewEdge(rSite a, rSite b) {
        voronoiDiagram.edges.emplace_back(make_unique<Edge>(a.point, b.point));
        return voronoiDiagram.edges.back().get();
    }

    Point* FortuneHyperbolicImplementation::getNewVertex(rPoint p) {
        voronoiDiagram.vertices.push_back(make_unique<Point>(p));
        return voronoiDiagram.vertices.back().get();
    }
}