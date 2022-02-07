#include <string>
#include "lib.hpp"

namespace hyperbolic {
    /*
    * writes a calculated Voronoi diagram and the corresponding sites to a file
     * */
    void draw_diagram(string& filename, VoronoiDiagram &voronoiDiagram, vector <Point> &sites);

    /*
     * writes the Delaunay triangulation of a Voronoi diagram to a file
     * */
    void write_delaunay(string& filename, VoronoiDiagram &voronoiDiagram);
}