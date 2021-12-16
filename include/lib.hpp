#pragma once
#include <vector>
#include <memory>

namespace hyperbolic {
    void print_version();

    typedef double float_t;

    struct Site {
        float_t theta, r;
        Site(float_t theta, float_t r) : theta(theta), r(r) {}
    };

    class VoronoiDiagram {
    public:
        virtual ~VoronoiDiagram() = default;
        virtual void calculate() = 0;
    };

    std::unique_ptr<VoronoiDiagram> getInstance(std::vector<Site> &sites);
}