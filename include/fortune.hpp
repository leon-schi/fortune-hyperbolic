#pragma once

#include <utility>
#include <vector>
#include <queue>
#include <memory>
#include <type_traits>

#include <lib.hpp>
#include <datastructures.hpp>
#include <beachline.hpp>
#include <kernels.hpp>

using std::vector, std::make_unique, std::make_shared;

namespace hyperbolic {
    template<class K>
    class FortuneHyperbolicImplementation {
    private:
            K kernel;
            BeachLine<K> beachLine;
            double r_sweep;

            EventQueue<SiteEvent> siteEventQueue;
            EventQueue<CircleEvent> circleEventQueue;

            vector<Site> sites;
            VoronoiDiagram& voronoiDiagram;

            bool isCalculated = false;

            bool eventsRemaining() {
                return !(siteEventQueue.empty() && circleEventQueue.empty());
            };

            static void invalidateCircleEvent(CircleEvent* e){
                if (e) e->invalidate();
            }

            void addCircleEvent(rBeachLineElement first, rBeachLineElement second) {
                invalidateCircleEvent(first.next.get());
                invalidateCircleEvent(second.previous.get());
                Point centerCircleEvent;
                if (kernel.predict_circle_event(centerCircleEvent, first.first, first.second, second.second)) {
                    _float_t radius = distance<_float_t>(centerCircleEvent, first.first.point);
                    shared_ptr<CircleEvent> ce = make_shared<CircleEvent>(first, second, centerCircleEvent, radius);
                    if (ce->r >= r_sweep) {
                        first.next = ce;
                        second.previous = ce;
                        circleEventQueue.push(ce);
                    }
                }
            };

            void initializeBeachLine() {
                // TODO: handle speicial case that the three closest sites are equidistant from the origin

                rSite a = siteEventQueue.top()->site;
                siteEventQueue.pop();
                rSite b = siteEventQueue.top()->site;
                siteEventQueue.pop();
                r_sweep = b.point.r;

                Edge* edge = getNewEdge(b, a, EdgeType::BIDIRECTIONAL);
                auto first = new BeachLineElement(b, a, edge, edge->firstVertex);
                auto last = new BeachLineElement(a, b, edge, edge->secondVertex);
                beachLine.insert(0, *first, *last);
            };

            void handleSiteEvent(const SiteEvent& e) {
#ifndef NDEBUG
                std::cout << "Handling Site Event with radius " << e.r << std::endl;
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
            };
            void handleCircleEvent(const CircleEvent& e) {
#ifndef NDEBUG
                if (!e.isValid())
                    std::cout << "Skipping invalid circle event" << std::endl;
                else
                    std::cout << "Handling Circle Event with radius " << e.r << " concerning beach line elements (" << e.first.first.ID << ", " << e.first.second.ID << ") and (" << e.second.first.ID << ", " << e.second.second.ID << ")" << std::endl;
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
            };

            Edge* getNewEdge(rSite a, rSite b, EdgeType edgeType) {
                voronoiDiagram.edges.emplace_back(make_unique<Edge>(a, b, edgeType));
                return voronoiDiagram.edges.back().get();
            };
            Point* getNewVertex(rPoint p) {
                voronoiDiagram.vertices.push_back(make_unique<Point>(p));
                return voronoiDiagram.vertices.back().get();
            };
        public:
            FortuneHyperbolicImplementation(VoronoiDiagram& v, const vector<Point>& sites) : voronoiDiagram(v), beachLine(kernel) {
                unsigned long long id = 0;
                for (const Point& p : sites) {
                    this->sites.emplace_back(Site(p, id));
                    id++;
                }
                r_sweep = 0;
            };

            void calculate() {
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
            };
    };
}