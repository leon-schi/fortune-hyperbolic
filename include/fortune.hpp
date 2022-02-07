#pragma once

#include <utility>
#include <vector>
#include <queue>
#include <memory>
#include <type_traits>

#include <lib.hpp>
#include <datastructures.hpp>

using std::vector;

namespace hyperbolic {

    class FortuneHyperbolicImplementation : public FortuneHyperbolic {
        private:
            BeachLine beachLine;
            _float_t r_sweep;

            EventQueue<SiteEvent> siteEventQueue;
            EventQueue<CircleEvent> circleEventQueue;
            bool eventsRemaining();

            VoronoiDiagram& voronoiDiagram;
            vector<Site> sites;

            bool isCalculated = false;

            void addCircleEvent(rBeachLineElement first, rBeachLineElement second);

            void initializeBeachLine();
            void sanitizeEdges();

            void handleSiteEvent(const SiteEvent& e);
            void handleCircleEvent(const CircleEvent& e);

            SiteTripleMap circleEventCache;
            bool predictCircleEvent(Point& result, rBeachLineElement a, rBeachLineElement b);

            Edge* getNewEdge(rSite a, rSite b, EdgeType edgeType);
            Point* getNewVertex(rPoint p);
        public:
            FortuneHyperbolicImplementation(VoronoiDiagram& v, const vector<Point>& sites);
            void calculate() override;
    };
}