import cairo
import math
import numpy as np

from utils import read_graph

line_width = 0.0005
point_radius = 0.0005
offset = 0.5
scale = 0.5
twopi = 2*math.pi

def read_sites(filename):
    coordinates = []
    with open(filename, 'r') as f:
        for line in f:
            theta, r = [float(x) for x in line.split(' ')]
            coordinates.append([r*math.cos(theta), r*math.sin(theta)])
    coordinates = np.array(coordinates)
    max_val = np.max(coordinates.ravel()) + 0.3
    coordinates /= max_val
    return coordinates

def transform(x, y):
    return x*scale+offset, y*scale+offset

def draw_sites(context, coordinates,):
    context.set_source_rgb(0, 0, 0)
    for x, y in coordinates:
        context.arc(*transform(x, y), point_radius, 0, twopi)
        context.fill()

def draw_edge(context, a, b, coordinates):
    context.move_to(*transform(*coordinates[a]))
    context.line_to(*transform(*coordinates[b]))
    context.stroke()

def draw_edges(context, adj1, adj2, coordinates):
    context.set_line_width(line_width)

    for u in adj1.keys():
        set1 = set(adj1[u])
        if u in adj2.keys():
            set2 = set(adj2[u])
        else:
            set2 = set()
        
        context.set_source_rgb(0.5, 0.5, 0.5)
        for v in set1.intersection(set2):
            draw_edge(context, u, v, coordinates)

        context.set_source_rgb(1, 0, 0)
        for v in set1 - set2:
            draw_edge(context, u, v, coordinates)

        context.set_source_rgb(0, 0, 1)
        for v in set2 - set1:
            draw_edge(context, u, v, coordinates)


def draw(coordinates, adj1, adj2):
    dim = 700
    with cairo.SVGSurface("out.svg", dim, dim) as surface:
        context = cairo.Context(surface)
        context.scale(dim, dim)

        draw_edges(context, adj1, adj2, coordinates)
        draw_sites(context, coordinates)


if __name__ == "__main__":
    adj1 = read_graph('../cmake-build-debug/bin/delaunay_lp.txt')
    adj2 = read_graph('../cmake-build-debug/bin/delaunay_hp.txt')
    coordinates = read_sites('../cmake-build-debug/bin/sample.txt')
    draw(coordinates, adj1, adj2)