#include <iostream>
#include <queue>

#include <lib.hpp>
#include <fortune.hpp>
#include <datastructures.hpp>
#include <cmake.hpp>

using std::cout, std::endl, std::unique_ptr, std::make_unique, std::vector;
using namespace hyperbolic;

namespace hyperbolic {
    void print_version() {
        cout << VERSION_MAJOR << endl;
    }

    unique_ptr<FortuneHyperbolic> getInstance(VoronoiDiagram& diagram, std::vector<Site>& sites) {
        return make_unique<FortuneHyperbolicImplementation>(diagram, sites);
    }

    void FortuneHyperbolicImplementation::handleSiteEvent(const SiteEvent * const e) {
        if (e == nullptr) return;

        cout << "Handling Site Event with radius " << e->getRadius() << endl;

        if (beachLine.size() == 0) {

        }
    }

    void FortuneHyperbolicImplementation::handleCircleEvent(const CircleEvent * const e) {
        if (e == nullptr) return;
        if (!e->isValid()) return;

        cout << "Handling Circle Event with radius " << e->getRadius() << endl;
    }

    bool FortuneHyperbolicImplementation::eventsRemaining() {
        return !(siteEventQueue.empty() && circleEventQueue.empty());
    }

    void FortuneHyperbolicImplementation::initializeBeachLine() {
        const SiteEvent* e1 = siteEventQueue.getTop();
        siteEventQueue.pop();
        const SiteEvent* e2 = siteEventQueue.getTop();
        siteEventQueue.pop();

        rSite a = e1->site, b = e2->site;
        Edge& edge = getNewEdge(a, b);
        auto first = new BeachLineElement(b, a, edge, edge.second);
        auto second = new BeachLineElement(a, b, edge, edge.first);
        beachLine.insert(0, *first, *second);
    }

    void FortuneHyperbolicImplementation::calculate() {
        if (isCalculated) return;
        isCalculated = true;

        if (sites.size() <= 1) return;

        for (rSite s : sites) {
            siteEventQueue.push(make_unique<SiteEvent>(s));
        }

        while (eventsRemaining()) {
            const SiteEvent* se = siteEventQueue.getTop();
            const CircleEvent* ce = circleEventQueue.getTop();

            if (!ce || (se && se->r <= ce->r)) {
                handleSiteEvent(se);
                siteEventQueue.pop();
            } else {
                handleCircleEvent(ce);
                circleEventQueue.pop();
            }
        }
    }

    Edge& FortuneHyperbolicImplementation::getNewEdge(rSite a, rSite b) {
        voronoiDiagram.edges.emplace_back(Edge(a, b));
        return voronoiDiagram.edges.back();
    }
}