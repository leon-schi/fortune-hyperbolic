#pragma once
#include <vector>
#include <memory>

#include <geometry.hpp>

using std::vector, std::unique_ptr;

namespace hyperbolic {
    void print_version();

    struct Edge {
        rPoint siteA, siteB;
        Point last_known_position;
        const Bisector bisector;
        _float_t theta_start = 0, theta_end = 0;
        pPoint firstVertex = nullptr, secondVertex = nullptr;
        Edge(rPoint a, rPoint b) : siteA(a), siteB(b), bisector(Bisector(&a, &b)) {};
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

    std::unique_ptr<FortuneHyperbolic> getNewInstance(VoronoiDiagram& diagram, const vector<Point>& sites);
    void draw_diagram(VoronoiDiagram& voronoiDiagram, vector<Point>& sites);
}