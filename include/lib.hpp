#pragma once
#include <vector>
#include <memory>

namespace hyperbolic {
    void print_version();

    typedef double _float_t;

#ifndef NDEBUG
    static int nextID = 0;
#endif

    struct Site {
#ifndef NDEBUG
        int ID;
#endif
        const _float_t theta, r;
        Site(_float_t theta, _float_t r) : theta(theta), r(r) {
#ifndef NDEBUG
            ID = nextID++;
#endif
        }
    };

    typedef const Site * pSite;
    typedef const Site & rSite;

    struct Vertex {
        _float_t r, theta;
        Vertex(_float_t r, _float_t theta) : r(r), theta(theta) {}
    };
    typedef const Vertex* pVertex;

    struct Edge {
        rSite a, b;
        pVertex first= nullptr, second = nullptr;
        Edge(rSite a, rSite b) : a(a), b(b) {};
    };

    class VoronoiDiagram {
    public:
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
    };

    class FortuneHyperbolic {
    public:
        virtual ~FortuneHyperbolic() = default;
        virtual void calculate() = 0;
    };

    std::unique_ptr<FortuneHyperbolic> getInstance(VoronoiDiagram& diagram, std::vector<Site>& sites);
}