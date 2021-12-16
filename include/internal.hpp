#pragma once

#include <utility>
#include <vector>
#include <queue>
#include "lib.hpp"

namespace hyperbolic {

    typedef const Site *const pSite;

    class Event {
        friend class EventQueue;
        protected:
            float_t r {0};
        public:
            virtual ~Event() = default;
            [[nodiscard]] float_t getRadius() const {
                return r;
            }
    };

    class SiteEvent : public Event {
    private:
        pSite site;
    public:
        explicit SiteEvent(pSite s) : Event(), site(s) {
            r = s->r;
        };
    };

    class CircleEvent : public Event {
        pSite a, b, c;
        void predict_radius();

    public:
        CircleEvent(pSite a, pSite b, pSite c) : Event(), a(a), b(b), c(c) {
            predict_radius();
        };
    };

    class EventQueue {
    private:
        struct CmpEvent {
            bool operator()(const Event* e1, const Event* e2) {
                return e1->r > e2->r;
            }
        };
        std::priority_queue<const Event* , std::vector<const Event*>, CmpEvent> q;
    public:
        void push(const Event* e) {
            q.push(e);
        }
        const Event* top() {
            return q.top();
        }
        void pop() {
            delete q.top();
            q.pop();
        }
        bool empty() {
            return q.empty();
        }
    };

    class FortuneHyperbolicImplementation : public VoronoiDiagram {
        private:
            std::vector<const Site *> sites;
            EventQueue q;

        public:
            explicit FortuneHyperbolicImplementation(std::vector<Site> &v) {
                for (auto &s : v)
                    sites.emplace_back(&s);
            };

            void calculate() override;
    };
}