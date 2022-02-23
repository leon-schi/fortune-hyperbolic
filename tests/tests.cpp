#include <gtest/gtest.h>

#include <fortune-hyperbolic/beachline.hpp>
#include <fortune-hyperbolic/kernels.hpp>
#include <fortune-hyperbolic/fortune.hpp>

#include <vector>
#include <memory>

using std::shared_ptr, std::make_shared;
using namespace hyperbolic;

TEST(BeachLineTest, InsertsCorrectly) {
    FullNativeKernel<double> K;
    BeachLine<FullNativeKernel<double>, double> beachLine(K);

    Point<double> mock(0, 0);
    Point<double>* pMock = &mock;

    int n = 30;
    vector<Site<double>> v;
    for (int i = 0; i < n; i++) {
        v.emplace_back(Point<double>(1, 1), i);
    }

    for (int i = 1; i < n; i++) {
        Site<double>* hitSite = &v[i];
        if (beachLine.size() > 0) {
            auto result = beachLine.getFirstElement();
            hitSite = &result->second;
        }
        auto* first = new BeachLineElement<double>(v[i], *hitSite, nullptr, pMock);
        auto* second = new BeachLineElement<double>(*hitSite, v[i], nullptr, pMock);
        beachLine.insert(0, *first, *second);
    }

    vector<BeachLineElement<double>*> elements;
    beachLine.getRemainingElements(elements);

    unsigned long long current_id = elements[0]->second.ID;
    for (size_t i = 1; i < elements.size(); i++) {
        EXPECT_EQ(current_id, elements[i]->first.ID);
        current_id = elements[i]->second.ID;
    }
};

TEST(VoronoiTest, ComputesCorrectly) {
    VoronoiDiagram v;
    vector<Point<double>> sites = {
            {3, 2.43},
            {2, 2.19},
            {6, 0.87},
            {9.2, 1.23}
    };
    FortuneHyperbolicImplementation<FullNativeKernel<double>, double>fortune(v, sites);
    fortune.calculate();

    EXPECT_EQ(3, v.edges.size());

    EXPECT_EQ(0, v.edges[0]->siteA.ID);
    EXPECT_EQ(1, v.edges[0]->siteB.ID);

    EXPECT_EQ(2, v.edges[1]->siteA.ID);
    EXPECT_EQ(1, v.edges[1]->siteB.ID);

    EXPECT_EQ(3, v.edges[2]->siteA.ID);
    EXPECT_EQ(1, v.edges[2]->siteB.ID);
}