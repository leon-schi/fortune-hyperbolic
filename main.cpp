#include <iostream>
#include <vector>
#include <lib.hpp>

using namespace std;
using namespace hyperbolic;

int main() {
    print_version();

    vector<Site> sites = {Site(4, 2.43), Site(2, 2.19), Site(7, 0.87)};
    auto v = getInstance(sites);

    v->calculate();
    return 0;
}