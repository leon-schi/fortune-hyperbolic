#include <iostream>
#include <queue>

#include <lib.hpp>
#include <internal.hpp>
#include <cmake.hpp>

using namespace std;
using namespace hyperbolic;

namespace hyperbolic {
    void print_version() {
        cout << VERSION_MAJOR << endl;
    }

    unique_ptr<VoronoiDiagram> getInstance(std::vector<Site> &sites) {
        return make_unique<FortuneHyperbolicImplementation>(sites);
    }

    void CircleEvent::predict_radius() {

    }

    void FortuneHyperbolicImplementation::calculate() {
        EventQueue q;
        for (const Site *s : sites) {
            q.push( new SiteEvent(s));
        }

        while (!q.empty()) {
            cout << q.top()->getRadius() << endl;
            q.pop();
        }
    }
}