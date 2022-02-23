#include <iostream>
#include <vector>
#include <chrono>

#include <fortune-hyperbolic/geometry.hpp>
#include <fortune-hyperbolic/canvas.hpp>
#include <fortune-hyperbolic/kernels.hpp>
#include <fortune-hyperbolic/fortune.hpp>

#include <boost/multiprecision/mpfr.hpp>

#include "cxxopts.hpp"

using namespace std;
using namespace hyperbolic;
using namespace boost;
using namespace multiprecision;

// use 200 bit floating point numbers
using floating_point_type = number<mpfr_float_backend<200, allocate_stack>>;
//using floating_point_type = double;

int main(int argc, char* argv[]) {

    cxxopts::Options options(
            argv[0], "Calculates a hyperbolic voronoi diagram from a set of sites using the hyperbolic version of Fortune's Algorithm. Version " + Utils::get_version());

    options.add_options()
            ("i,input", "Input Filename", cxxopts::value<std::string>())
            ("v,verbose", "Enable verbose output (only for debugging)", cxxopts::value<bool>()->default_value("false"))
            ("d,output_diagram", "Output Filename for writing the diagram svg", cxxopts::value<std::string>())
            ("t,output_triangulation", "Output Filename for writing the delaunay triangulation", cxxopts::value<std::string>())
            ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    string input_file = result["i"].as<string>();

    // read the input
    ifstream input_stream(input_file);
    vector<Point<double>> sites;
    try {
        if (input_stream.is_open()) {
            string line;
            while (getline(input_stream, line)) {
                _float_t theta, r;
                auto pos = line.find(' ');
                theta = stod(line.substr(0, pos));
                r = stod(line.substr(pos));
                sites.emplace_back(r, theta);
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

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    VoronoiDiagram v;
    FortuneHyperbolicImplementation<FullNativeKernel<floating_point_type>, floating_point_type> fortune(v, sites, result["v"].as<bool>());
    fortune.calculate();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    cout << "Finished calculating Voronoi diagram after " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    VoronoiCanvasOptions canvas_options;
    canvas_options.width = 500;
    VoronoiCanvas canvas(v, sites);
    if (result.count("d")) {
        string output_file = result["d"].as<string>();
        canvas.set_options(canvas_options);
        canvas.draw_diagram(output_file);
    }

    if (result.count("t")) {
        string output_file = result["t"].as<string>();
        canvas.write_delaunay_triangulation(output_file);
    }
    return 0;
}