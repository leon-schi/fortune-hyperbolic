#pragma once
#include <vector>
#include <memory>

namespace hyperbolic {
    void print_version();

    typedef double _float_t;

    struct Point {
        _float_t r = 0, theta = 0;
        Point() = default;
        Point(_float_t r, _float_t theta) : r(r), theta(theta) {}
    };
    typedef const Point* pPoint;
    typedef const Point& rPoint;

    struct Edge {
        rPoint siteA, siteB;
        _float_t theta_start = 0, theta_end = 0;
        pPoint firstVertex= nullptr, secondVertex = nullptr;
        Edge(rPoint a, rPoint b) : siteA(a), siteB(b) {};
    };

    class VoronoiDiagram {
    public:
        std::vector<Point> vertices;
        std::vector<Edge> edges;
    };

    class FortuneHyperbolic {
    public:
        virtual ~FortuneHyperbolic() = default;
        virtual void calculate() = 0;
    };

    std::unique_ptr<FortuneHyperbolic> getInstance(VoronoiDiagram& diagram, const std::vector<Point>& sites);
}