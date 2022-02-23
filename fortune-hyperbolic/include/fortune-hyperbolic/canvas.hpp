#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <list>
#include <cmath>
#include <vector>

#include <fortune-hyperbolic/geometry.hpp>
#include <fortune-hyperbolic/calculations.hpp>

using std::string, std::max_element, std::list;

namespace hyperbolic {
    /**
    * Drawing options that can be provided to the VoronoiCanvas class
    * */
    struct VoronoiCanvasOptions {
        double width = 300, height = 300;
        double line_width = 0.01, point_width = 0.02;
        double resolution = 0.01;
        string voronoi_edge_color = "black";
        string delaunay_edge_color = "red";
        string point_color = "black";
        bool draw_delaunay = true;
    };

    struct CartesianPoint {
        double x, y;
        CartesianPoint() = default;
        explicit CartesianPoint(Point<double> p) {
            x = cos(p.theta)*p.r;
            y = sin(p.theta)*p.r;
        }
        CartesianPoint(const CartesianPoint& p) = default;
        CartesianPoint(double x, double y) : x(x), y(y) {}
        CartesianPoint operator *(const double scale) const {
            return {x*scale, y*scale};
        }
    };

    using Path = vector<CartesianPoint>;

    /**
     * Class used for drawing Voronoi diagrams and Delaunay triangulations
     * */
    class VoronoiCanvas {
    private:
        VoronoiDiagram& voronoiDiagram;
        vector<Point<double>>& sites;

        VoronoiCanvasOptions options;

        double max_r = 0, scale = 0;
        CartesianPoint offset = {0, 0};

        static void add_delaunay_edge(const Edge& e, Path& p) {
            p.push_back(CartesianPoint(e.siteA.point));
            p.push_back(CartesianPoint(e.siteB.point));
        }

        bool out_of_bounds(Point<double>& p) const {
            return
                (std::abs(cos(p.theta)*p.r*scale) + offset.x > options.width) ||
                (std::abs(sin(p.theta)*p.r*scale) + offset.y > options.height);
        }

        void render_edge(Point<double>* from, Point<double>* to, HyperboloidBisector<double>& b, list<CartesianPoint>& p, bool ccw=true) const {
            double dt = options.resolution * scale;
            Point u_native(b.u);
            double t = (from) ? distance<double>(*from, u_native) : 0;
            double t_end = (to) ? distance<double>(*to, u_native) : DBL_MAX;

            HyperboloidVec<double> v = b.v;
            Point v_polar(v);
            Point u_polar(b.u);
            double theta = clip(v_polar.theta - u_polar.theta);
            if ((ccw && theta >= M_PI) || (!ccw && theta <= M_PI))
                v = v*(-1);

            while (t < t_end) {
                Point point(b.u*cosh(t) + v*sinh(t));
                p.emplace_back(point);
                t += dt;
                if (out_of_bounds(point)) break;
                //if (point.r*scale*1.414 > width) break;
            }

            if (from)
                p.emplace_front(*from);
            if (to)
                p.emplace_back(*to);
        }

        void add_edge(Edge& e, Path& p) const {
            HyperboloidBisector<double> b(e.siteA.point, e.siteB.point);
            Point q(b.u);
            if (e.edgeType == EdgeType::BIDIRECTIONAL) {
                list<CartesianPoint> cw, ccw;
                render_edge(nullptr, e.firstVertex, b, ccw, true);
                render_edge(nullptr, e.secondVertex, b, cw, false);
                cw.emplace_front(q);
                cw.reverse();
                p.insert(p.begin(), ccw.begin(), ccw.end());
                p.insert(p.begin(), cw.begin(), cw.end());
            }
            else {
                list<CartesianPoint> p_list;
                render_edge(e.firstVertex, e.secondVertex, b, p_list, e.edgeType == EdgeType::CCW);
                p.insert(p.begin(), p_list.begin(), p_list.end());
            }
        }

        void save_to_file(string& filename) const {
            string representation;
            to_svg_representation(representation);

            std::fstream output_file_stream(filename, std::fstream::out);
            output_file_stream << representation;
        }

        void to_svg_representation(std::string &svg_representation) const {

            svg_representation =
                    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE svg PUBLIC "
                    "\"-//W3C//DTD SVG 1.1//EN\" "
                    "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n<svg "
                    "xmlns=\"http://www.w3.org/2000/svg\"\nxmlns:xlink=\"http://www.w3.org/"
                    "1999/xlink\" "
                    "xmlns:ev=\"http://www.w3.org/2001/xml-events\"\nversion=\"1.1\" ";

            svg_representation += std::string("baseProfile=\"full\"\nwidth=\"") +
                                  std::to_string(options.width) + "\" height=\"" +
                                  std::to_string(options.height) + "\">\n\n";


            for (auto& e : voronoiDiagram.edges) {
                std::string path_representation;
                Path p;

                add_edge(*e, p);
                svg_path_representation(p, path_representation, offset, scale, options.voronoi_edge_color);
                svg_representation += path_representation;

                if (options.draw_delaunay) {
                    p.clear();
                    add_delaunay_edge(*e, p);
                    svg_path_representation(p, path_representation, offset, scale, options.delaunay_edge_color);
                    svg_representation += path_representation;
                }
            }

            for (const Point<double>& s : sites) {
                std::string point_representation;
                svg_point_representation(CartesianPoint(s), options.point_width, point_representation, options.point_color);
                svg_representation += point_representation;
            }

            std::string point_representation;
            svg_point_representation(CartesianPoint(0, 0), options.point_width/2, point_representation, options.point_color);
            svg_representation += point_representation;

            svg_representation += "\n</svg>\n";
        }

        static void svg_path_representation(const Path& path, string& path_representation, const CartesianPoint& offset, double scale, const string& color="black") {
            path_representation = string("");

            if (!path.empty()) {
                path_representation = std::string("<path d =\"");
                // Print the first point of the path.
                CartesianPoint point(path.front()*scale);
                point.x += offset.x;
                point.y += offset.y;
                path_representation += std::string("M ") + std::to_string(point.x) +
                                       "," + std::to_string(point.y) + " ";

                //Print the remaining points.
                for (size_t index = 1; index < path.size(); ++index) {
                    point = CartesianPoint(path[index]*scale);
                    point.x += offset.x;
                    point.y += offset.y;
                    path_representation += std::string("L ") + std::to_string(point.x) + ", " + std::to_string(point.y) + " ";
                }

                double path_width = 0.01 * scale;
                path_representation +=
                        std::string("\" stroke = \"") + color + "\" stroke-width = \"" +
                        std::to_string(path_width) + R"(" fill="none"/>)";
            }
        }

        void svg_point_representation(const CartesianPoint& p, double radius, std::string &svg_point_representation, const string& color) const {
            CartesianPoint center(p*scale);
            center.x += offset.x;
            center.y += offset.y;

            double stroke_width = options.line_width * scale;
            svg_point_representation =
                    std::string("<circle cx=\"") + std::to_string(center.x) + "\" cy=\"" +
                    std::to_string(center.y) + "\" r=\"" +
                    std::to_string(radius * scale) + "\" fill=\"" + color +
                    "\" stroke=\"" + color + "\" stroke-width=\"" +
                    std::to_string(stroke_width) + "\"/>\n";
        }

    public:
        explicit VoronoiCanvas(VoronoiDiagram& v, vector<Point<double>>& sites) : voronoiDiagram(v), sites(sites) {};

        void set_options(VoronoiCanvasOptions& opt) {
            options = opt;
        }

        /**
        * writes a calculated Voronoi diagram and the corresponding sites to a file
        * */
        void draw_diagram(string filename) {
            max_r = (max_element(sites.begin(), sites.end()))->r;
            double min_hw = min(options.width, options.height);
            scale = min_hw/(2*max_r);
            offset = CartesianPoint(options.width/2, options.height / 2);
            save_to_file(filename);
        }

        /**
         * writes the Delaunay triangulation of a Voronoi diagram to a file
         * */
        void write_delaunay_triangulation(string filename) {
            std::fstream output_file_stream(filename, std::fstream::out);
            for (auto& e : voronoiDiagram.edges) {
                output_file_stream << e->siteA.ID << " " << e->siteB.ID << "\n";
            }
        }
    };
}