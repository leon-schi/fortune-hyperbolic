#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>

#include <lib.hpp>
#include <geometry.hpp>
#include <utils.hpp>

#include <kernels.hpp>

#include "cxxopts.hpp"

using namespace std;
using namespace hyperbolic;

int main(int argc, char* argv[]) {

    cxxopts::Options options(
            argv[0], "Calculates a hyperbolic voronoi diagram from a set of sites using the hyperbolic version of Fortune's Algorithm. Version " + get_version());

    options.add_options()
            ("i,input", "Input Filename", cxxopts::value<std::string>())
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
    vector<Point> sites;
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
    auto fortune = getNewInstance(v, sites);
    fortune->calculate();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    cout << "Finished calculating Voronoi diagram after " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    if (result.count("d")) {
        string output_file = result["d"].as<string>();
        draw_diagram(output_file, v, sites);
    }

    if (result.count("t")) {
        string output_file = result["t"].as<string>();
        write_delaunay(output_file, v);
    }
    return 0;
}