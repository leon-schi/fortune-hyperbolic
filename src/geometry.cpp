#include <cmath>
#include <vector>

#include <lib.hpp>
#include <calculations.hpp>
#include <iostream>
#include <fstream>

using std::vector;

namespace hyperbolic {
    struct CartesianPoint {
        double x, y;
        explicit CartesianPoint(Point p) {
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

    class VoronoiCanvas {
    public:
        VoronoiDiagram& voronoiDiagram;
        vector<Point>& sites;
        double width = 300;

        explicit VoronoiCanvas(VoronoiDiagram& v, vector<Point>& sites) : voronoiDiagram(v), sites(sites) {};

        /*
         * Transforms the bisector b in the interval [start, end] to a Path.
         *  resolution is an angle that describes the steps in angular coordinates at which points are sampled from the bisector
         * */
        void add_path_from_straight_edge(const Edge& e, Path& p) const {
            const Bisector& b = e.bisector;
            if (e.firstVertex && e.secondVertex) {
                p.push_back(CartesianPoint(*e.firstVertex));
                p.push_back(CartesianPoint(*e.secondVertex));
            } else if (e.firstVertex || e.secondVertex) {
                pPoint v = (e.firstVertex) ? e.firstVertex : e.secondVertex;
                p.push_back(CartesianPoint(*v));
                p.push_back(CartesianPoint(Point(width*width, e.last_known_position.theta)));
            } else {
                p.push_back(CartesianPoint(Point(width*width, b.straight_angle + M_PI)));
                p.push_back(CartesianPoint(Point(width*width, b.straight_angle)));
            }
        }

        void add_path_from_curved_edge(const Edge& e, double resolution, Path& p) const {
            const Bisector& b = e.bisector;
            double start = b.theta_start, end = b.theta_end;

            if (e.firstVertex && e.secondVertex) {
                start = e.firstVertex->theta; end = e.secondVertex->theta;
                if (clip(start - b.theta_start) > clip(end - b.theta_start))
                    std::swap(start, end);
            } else if (e.firstVertex || e.secondVertex) {
                pPoint v = (e.firstVertex) ? e.firstVertex : e.secondVertex;
                if (clip(v->theta - b.theta_start) < clip(e.last_known_position.theta - b.theta_start)) {
                    start = v->theta;
                    end = b.theta_end-resolution;
                } else {
                    start = b.theta_start+resolution;
                    end = v->theta;
                }
            } else {
                start = b.theta_start+resolution; end = b.theta_end-resolution;
            }

            end = clip(end - start);
            double current = 0;
            while (current <= end) {
                p.push_back(CartesianPoint(Point(b(current + start), current + start)));
                current += resolution;
            }
    }

        void save_to_file(string& filename) const {
            string representation;
            to_svg_representation(representation);

            std::fstream output_file_stream(filename, std::fstream::out);
            output_file_stream << representation;
        }

        void to_svg_representation(std::string &svg_representation) const {

            //TODO: find the site with the largest radius to obtain good bounds

            double scale = width/(2*3);
            CartesianPoint offset(width/2, width/2);

            svg_representation =
                    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE svg PUBLIC "
                    "\"-//W3C//DTD SVG 1.1//EN\" "
                    "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n<svg "
                    "xmlns=\"http://www.w3.org/2000/svg\"\nxmlns:xlink=\"http://www.w3.org/"
                    "1999/xlink\" "
                    "xmlns:ev=\"http://www.w3.org/2001/xml-events\"\nversion=\"1.1\" ";

            svg_representation += std::string("baseProfile=\"full\"\nwidth=\"") +
                                  std::to_string(width) + "\" height=\"" +
                                  std::to_string(width) + "\">\n\n";


            for (auto& e : voronoiDiagram.edges) {
                std::string path_representation;
                Path p;
                if (e->bisector.is_straight)
                    add_path_from_straight_edge(*e, p);
                else
                    add_path_from_curved_edge(*e, 0.0001, p);
                svg_path_representation(p, path_representation, offset, scale);
                svg_representation += path_representation;
            }

            for (const Point& s : sites) {
                std::string point_representation;
                svg_point_representation(CartesianPoint(s), 0.02, point_representation, offset, scale);
                svg_representation += point_representation;
            }

            std::string point_representation;
            svg_point_representation(CartesianPoint(0, 0), 0.01, point_representation, offset, scale);
            svg_representation += point_representation;

            svg_representation += "\n</svg>\n";
        }

        static void svg_path_representation(const Path& path, string& path_representation, const CartesianPoint& offset, double scale) {
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

                double path_width = 0.02 * scale;
                path_representation +=
                        std::string("\" stroke = \"") + "black" + "\" stroke-width = \"" +
                        std::to_string(path_width) + R"(" fill="none"/>)";
            }
        }

        static void svg_point_representation(const CartesianPoint& p, double radius, std::string &svg_point_representation, const CartesianPoint& offset, double scale) {
            CartesianPoint center(p*scale);
            center.x += offset.x;
            center.y += offset.y;

            double stroke_width = 0.02 * scale;
            svg_point_representation =
                    std::string("<circle cx=\"") + std::to_string(center.x) + "\" cy=\"" +
                    std::to_string(center.y) + "\" r=\"" +
                    std::to_string(radius * scale) + "\" fill=\"" + "black" +
                    "\" stroke=\"" + "black" + "\" stroke-width=\"" +
                    std::to_string(stroke_width) + "\"/>\n";
        }
    };

    void draw_diagram(VoronoiDiagram& voronoiDiagram, vector<Point>& sites) {
        VoronoiCanvas canvas(voronoiDiagram, sites);
        string filename = "diagram.svg";
        canvas.save_to_file(filename);
    };
}