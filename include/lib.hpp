#pragma once
#include <vector>
#include <memory>
#include <string>

#include <geometry.hpp>

using std::vector, std::unique_ptr, std::string;

namespace hyperbolic {
    string get_version();

    enum EdgeType {CCW, CW, BIDIRECTIONAL};

    /*
     * represents an edge of the Voronoi diagram
     * */
    struct Edge {
        EdgeType edgeType;
        // the sites that define the bisector corresponding to the edge
        rSite siteA, siteB;
        // A point that is on the edge
        Point last_known_position;
        // the bisector corresponding to the edge
        const Bisector bisector;
        // the vertices of the Voronoi diagram the edge is adjacent to
        pPoint firstVertex = nullptr, secondVertex = nullptr;
        Edge(rSite a, rSite b, EdgeType edgeType) :
            edgeType(edgeType), siteA(a), siteB(b), bisector(Bisector(&a.point, &b.point)) {};
    };

    class VoronoiDiagram {
    public:
        vector<unique_ptr<Point>> vertices;
        vector<unique_ptr<Edge>> edges;
    };

    class FortuneHyperbolic {
    public:
        virtual ~FortuneHyperbolic() = default;
        virtual void calculate() = 0;
    };


    /*
     * returns an instance of the class that implements the algorithm
     * */
    std::unique_ptr<FortuneHyperbolic> getNewInstance(VoronoiDiagram& diagram, const vector<Point>& sites);
}