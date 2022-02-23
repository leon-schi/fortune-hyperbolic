#pragma once

#include <cmath>
#include <vector>
#include <memory>

using std::vector, std::max, std::min, std::unique_ptr;

namespace hyperbolic {
    typedef double _float_t;

    template<typename _float_T>  struct Point;
    template<typename _float_T> struct HyperboloidVec;

    /*
     * Most basic class for representing a point in polar coordinates
     * */
    template<typename _float_T>
    struct Point {
        _float_T r = 0, theta = 0;
        Point() = default;
        Point(_float_T r, _float_T theta) : r(r), theta(theta) {}
        explicit Point(const HyperboloidVec<_float_T>& p);
        bool operator < (const Point<_float_T>& p) const {
            return r < p.r;
        }
        explicit operator Point<double> () const {
            return Point<double>(static_cast<double>(r), static_cast<double>(theta));
        }
    };
    /*
     * Sites for internal use that assign an ID to a point which is used for caching
     * */
    template<typename _float_T>
    struct Site {
        Point<_float_T> point;
        const unsigned long long ID;
        Site(Point<_float_T> point, unsigned long long id) : point(point), ID(id) {};
        explicit operator Site<double> () const {
            return Site<double>(static_cast<Point<double>>(point), ID);
        }
    };

    enum EdgeType {CCW, CW, BIDIRECTIONAL};

    /*
     * represents an edge of the Voronoi diagram
     * */
    struct Edge {
        EdgeType edgeType;
        // the sites that define the bisector corresponding to the edge
        Site<double> siteA, siteB;
        // the vertices of the Voronoi diagram the edge is adjacent to
        Point<double> *firstVertex = nullptr, *secondVertex = nullptr;
        Edge(Site<double> a, Site<double> b, EdgeType edgeType) :
                edgeType(edgeType), siteA(a), siteB(b) {};
    };

    class VoronoiDiagram {
    public:
        vector<unique_ptr<Point<double>>> vertices;
        vector<unique_ptr<Edge>> edges;
    };
}