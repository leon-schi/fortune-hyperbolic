#include <iostream>
#include <vector>
#include <lib.hpp>
#include <geometry.hpp>

using namespace std;
using namespace hyperbolic;

int main() {
    print_version();

    vector<Point> sites = {
            Point(2, 0),
            Point(1.5, 1),
            Point(2.1, 0.6),
            Point(2, 3),
            Point(2, 5),
            Point(2.5, 3.5),
            Point(2.5, 2.2),
            Point(2, 2.5),
            Point(2, 1.5),
            Point(2, 1.7),
            Point(1, 5)
    };
    VoronoiDiagram v;
    auto fortune = getNewInstance(v, sites);

        fortune->calculate();

    draw_diagram(v, sites);

    return 0;
}