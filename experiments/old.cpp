#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <string>
#include <fstream>
#include <random>

#include <cxxopts.hpp>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Gt;
typedef Gt::Point_2                                         Point_2;
typedef Gt::FT                                              FT;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt>       Dt;
typedef CGAL::Creator_uniform_2<FT, Point_2>                Creator;

using namespace std;

int main(int argc, char** argv) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    CGAL::Timer timer;
    CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc;
    std::vector<Point_2> pts;

    cxxopts::Options options(
            argv[0], "Calculates a hyperbolic voronoi diagram from a set of sites using the implementation of CGAL 5.2 operating in the Poincare disk model");

    options.add_options()
            ("i,input", "Input Filename", cxxopts::value<std::string>())
            ("o,output", "Output Filename", cxxopts::value<std::string>())
            ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    string input_file = result["i"].as<string>();
    string output_file = result["o"].as<string>();

    // read the input
    ifstream input_stream(input_file);
    try {
        if (input_stream.is_open()) {
            string line;
            while (getline(input_stream, line)) {
                double theta, r;
                auto pos = line.find(' ');
                theta = stod(line.substr(0, pos));
                r = stod(line.substr(pos));
                r = tanh(r/2);
                Point_2 p(r*cos(theta), r*sin(theta));
                pts.push_back(p);
            }
            input_stream.close();
        } else {
            cout << "Unable to open input file \"" << input_file << "\"\n";
            exit(1);
        }
    } catch (...) {
        cout << "Error while reading file \"" << input_file << "\"\n";
        exit(1);
    }

    unordered_map<Dt::Vertex_handle, size_t> vertext_to_id;
    Dt dt_end;
    std::cout << "Insertion of point set (hyperbolic filtering only once at the end)"  << std::endl;
    std::cout << "===================================================================" << std::endl;
    timer.reset();
    timer.start();
    for (size_t i = 0; i < pts.size(); i++) {
        Dt::Vertex_handle h = dt_end.insert(pts[i]);
        vertext_to_id[h] = i;
    }
    timer.stop();
    std::cout << "Number of vertices:         " << dt_end.number_of_vertices() << std::endl;
    std::cout << "Number of hyperbolic faces: " << dt_end.number_of_hyperbolic_faces() << std::endl;
    std::cout << "Number of hyperbolic edges: " << dt_end.number_of_hyperbolic_edges() << std::endl;
    std::cout << "Time:                       " << timer.time() << std::endl;


    std::fstream output_file_stream(output_file, std::fstream::out);
    for (Dt::All_edges_iterator e = dt_end.all_edges_begin(); e != dt_end.all_edges_end(); e++) {
        Dt::Face_handle f = e->first;
        Dt::Vertex_handle a = f->vertex(f->cw(e->second));
        Dt::Vertex_handle b = f->vertex(f->ccw(e->second));
        output_file_stream << vertext_to_id[a->handle()] << " " << vertext_to_id[b->handle()] << "\n";
    }

    return EXIT_SUCCESS;
}
