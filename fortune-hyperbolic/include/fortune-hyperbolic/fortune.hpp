#pragma once

#include <utility>
#include <vector>
#include <queue>
#include <memory>
#include <type_traits>

#include <fortune-hyperbolic/geometry.hpp>
#include <fortune-hyperbolic/datastructures.hpp>
#include <fortune-hyperbolic/beachline.hpp>
#include <fortune-hyperbolic/kernels.hpp>
#include <cmake.hpp>

using std::vector, std::make_unique, std::make_shared, std::to_string;

namespace hyperbolic {
    class Utils {
    public:
        static string get_version() {
            return to_string(VERSION_MAJOR) + "." + to_string(VERSION_MINOR);
        }
    };

    /**
     * Main class that implements the algorithm. Requires a kernel K and a floating point type _float_T to use.
     * The kernel floating point type must match _float_T
     * */
    template<class K, typename _float_T>
    class FortuneHyperbolicImplementation {
    private:
            bool verbose;

            K kernel;
            BeachLine<K, _float_T> beachLine;
            _float_T r_sweep;

            using pBeachLineElement = BeachLineElement<_float_T>*;
            using rBeachLineElement = BeachLineElement<_float_T>&;
            using rSite = Site<_float_T>&;
            using pSite = Site<_float_T>*;

            EventQueue<SiteEvent<_float_T>, _float_T> siteEventQueue;
            EventQueue<CircleEvent<_float_T>, _float_T> circleEventQueue;

            vector<Site<_float_T>> sites;
            VoronoiDiagram& voronoiDiagram;

            bool isCalculated = false;

            bool eventsRemaining() {
                return !(siteEventQueue.empty() && circleEventQueue.empty());
            };

            static void invalidateCircleEvent(CircleEvent<_float_T>* e){
                if (e) e->invalidate();
            }

            /**
             * method that predicts and assigns new circle events
             * */
            void addCircleEvent(rBeachLineElement first, rBeachLineElement second) {
                invalidateCircleEvent(first.next.get());
                invalidateCircleEvent(second.previous.get());
                Point<_float_T> centerCircleEvent;
                if (kernel.predict_circle_event(centerCircleEvent, first.first, first.second, second.second)) {
                    _float_T radius = distance<_float_T>(centerCircleEvent, first.first.point);
                    shared_ptr<CircleEvent<_float_T>> ce = make_shared<CircleEvent<_float_T>>(first, second, centerCircleEvent, radius);
                    if (ce->r >= r_sweep) {
                        first.next = ce;
                        second.previous = ce;
                        circleEventQueue.push(ce);
                    }
                }
            };

            /**
             * Initializes the first elements in the beach line.
             * */
            void initializeBeachLine() {
                rSite a = siteEventQueue.top()->site;
                siteEventQueue.pop();
                rSite b = siteEventQueue.top()->site;
                siteEventQueue.pop();
                r_sweep = b.point.r;

                Edge* edge = getNewEdge(b, a, EdgeType::BIDIRECTIONAL);
                auto first = new BeachLineElement<_float_T>(b, a, edge, edge->firstVertex);
                auto last = new BeachLineElement<_float_T>(a, b, edge, edge->secondVertex);
                beachLine.insert(0, *first, *last);
            };

            /**
             * Called when a site event is handled.
             * */
            void handleSiteEvent(const SiteEvent<_float_T>& e) {
#ifndef NDEBUG
                if (verbose)
                    std::cout << "Handling Site Event with radius " << e.r << std::endl;
#endif
                r_sweep = e.r;

                pBeachLineElement first, second;
                int positionFirst = beachLine.search(e.site.point, e.r, first, second);
                rSite hitSite = first->second;

                Edge* edge = getNewEdge(e.site, hitSite, EdgeType::BIDIRECTIONAL);
                auto firstNew = new BeachLineElement<_float_T>(e.site, hitSite, edge, edge->firstVertex);
                auto secondNew = new BeachLineElement<_float_T>(hitSite, e.site, edge, edge->secondVertex);

                addCircleEvent(*first, *secondNew);
                addCircleEvent(*firstNew, *second);

                beachLine.insert(positionFirst, *firstNew, *secondNew);
            };

            /**
             * Called when a circle event is handled.
             * */
            void handleCircleEvent(const CircleEvent<_float_T>& e) {
#ifndef NDEBUG
                if (verbose) {
                    if (!e.isValid())
                        std::cout << "Skipping invalid circle event" << std::endl;
                    else
                        std::cout << "Handling Circle Event with radius " << e.r << " concerning beach line elements ("
                                  << e.first.first.ID << ", " << e.first.second.ID << ") and (" << e.second.first.ID
                                  << ", " << e.second.second.ID << ")" << std::endl;
                }
#endif

                if (!e.isValid()) return;
                r_sweep = e.r;
                Point<double>* v = getNewVertex(e.center);

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
                voronoiDiagram.edges.emplace_back(make_unique<Edge>(static_cast<Site<double>>(a), static_cast<Site<double>>(b), edgeType));
                return voronoiDiagram.edges.back().get();
            };
            Point<double>* getNewVertex(const Point<_float_T>& p) {
                voronoiDiagram.vertices.push_back(make_unique<Point<double>>(p));
                return voronoiDiagram.vertices.back().get();
            };
        public:
            /**
             * Instantiate the class.
             *
             * @param v Reference to the Voronoi Diagram to be computed
             * @param sites A vector of sites used as input
             * @param verbose Whether to print the beach line in each iteration or not
             */
            FortuneHyperbolicImplementation(VoronoiDiagram& v, const vector<Point<double>>& sites, bool verbose=false) : verbose(verbose), beachLine(kernel), voronoiDiagram(v) {
                unsigned long long id = 0;
                for (const Point<double>& p : sites) {
                    this->sites.emplace_back(Site<_float_T>(Point<_float_T>(p.r, p.theta), id));
                    id++;
                }
                r_sweep = 0;
            };

            /**
             * calculates the Voronoi diagram.
             * */
            void calculate() {
                if (isCalculated) return;
                isCalculated = true;

                if (sites.size() <= 1) return;

                for (rSite s : sites) {
                    siteEventQueue.push(make_shared<SiteEvent<_float_T>>(s));
                }

                initializeBeachLine();
#ifndef NDEBUG
                if (verbose) beachLine.print(r_sweep);
#endif

                while (eventsRemaining()) {
                    const SiteEvent<_float_T>* se = siteEventQueue.getTop();
                    const CircleEvent<_float_T>* ce = circleEventQueue.getTop();

                    if (!ce || (se && CmpEvent<Event<_float_T>, _float_T>::greater_equal(ce, se))) {
                        handleSiteEvent(*se);
                        siteEventQueue.pop();
                    } else {
                        handleCircleEvent(*ce);
                        circleEventQueue.pop();
                    }

#ifndef NDEBUG
                    if (verbose) beachLine.print(r_sweep);
#endif
                }
            };
    };
}