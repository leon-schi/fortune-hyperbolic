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

            EventQueue<SiteEvent> siteEventQueue;
            EventQueue<CircleEvent> circleEventQueue;

            VoronoiDiagram& voronoiDiagram;
            vector<Site> sites;

            bool isCalculated = false;

            void addCircleEvent(rBeachLineElement first, rBeachLineElement second);

            void handleSiteEvent(const SiteEvent& e);
            void handleCircleEvent(const CircleEvent& e);
            void initializeBeachLine();

            bool eventsRemaining();
            Edge& getNewEdge(rSite a, rSite b);
            Point& getNewVertex(rPoint p);
        public:
            FortuneHyperbolicImplementation(VoronoiDiagram& v, const vector<Point>& sites);

            void calculate() override;
    };
}