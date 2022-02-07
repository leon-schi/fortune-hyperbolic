#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>

#include "cxxopts.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    cxxopts::Options options(
            "generator", "Generator for sampling points in the hyperbolic plane and writing them to a file.");

    options.add_options()
            ("o,output", "Output Filename", cxxopts::value<std::string>()->default_value("sample.txt"))
            ("N", "Number of points to sample", cxxopts::value<int>()->default_value("100")) // a bool parameter
            ("R", "Radius within which points are sampled", cxxopts::value<double>()->default_value("5"))
            ("a,alpha", "Parameter alpha of the distribution from which we sample",cxxopts::value<double>()->default_value("1"))
            ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    int N = result["N"].as<int>();
    string filename = result["o"].as<string>();
    double R = result["R"].as<double>();
    double alpha = result["alpha"].as<double>();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::fstream output_file_stream(filename, std::fstream::out);
    output_file_stream << std::fixed << std::setprecision(6);

    for(int i=0; i<N; ++i) {
        double angle = dis(gen) * 2 * M_PI;
        double radius = acosh(1+(cosh(alpha*R)-1)*dis(gen))/alpha;

        output_file_stream << angle << " " << radius << "\n";
    }
    return 0;
}