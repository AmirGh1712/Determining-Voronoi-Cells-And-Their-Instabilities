import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
import random
import json
import sys

class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # creates a point from a given verctor
    def from_vector(v):
        return Point(v.x, v.y, v.z)

    # adds coordinates of current point with the coordinates of input point
    def __add__(self, o):
        return Point(self.x + o.x, self.y + o.y, self.z + o.z)

    # subtract the coordinates of the current point by the coordinates of the input point
    def __sub__(self, o):
        return Point(self.x - o.x, self.y - o.y, self.z - o.z)

    # multiply coordinates of current point with the coordinates of input point
    def mul(self, o):
        if (type(o) == int):
            return Point(self.x * o, self.y * o, self.z * o)

    # divide the coordinates of the current point by the coordinates of the input point
    def div(self, o):
        if (type(o) == int):
            return Point(self.x / o, self.y / o, self.z / o)

    # return true if the points have identical coordinates
    def __eq__(self, o):
        if o == None:
            return False
        # if not issubclass(type(o), Point):
        #     print("not type")
        #     return False
        return self.x == o.x and self.y == o.y and self.z == o.z

    # returns a string the represents the point
    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"



class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # adds coordinates of the current vector with the coordinates of input vector
    def __add__(self, o):
        return Vector(self.x + o.x, self.y + o.y, self.z + o.z)

    # multiply coordinates of the current vector with the coordinates of input vector
    def mul(self, o):
        if (type(o) == int):
            return Vector(self.x * o, self.y * o, self.z * o)

    # return true if the vectors have identical coordinates
    def __eq__(self, o):
        return round(self.x-o.x,4) == 0 and round(self.y-o.y,4) == 0 and round(self.z-o.z,4)==0

    # returns a string the represents the vector
    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z)  + ")->"


class Plane:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    # calculates a plane by given three points (in 3D)
    def three_points_to_plane(p1, p2, p3):
        a = (p2.y - p1.y) * (p3.z - p1.z) - (p3.y - p1.y) * (p2.z - p1.z)
        b = (p2.z - p1.z) * (p3.x - p1.x) - (p3.z - p1.z) * (p2.x - p1.x)
        c = (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)
        d = -p1.x * a - p1.y * b - p1.z * c
        return Plane(a, b, c, d)

    # calculates a plane by a given point in the vector and the normal vector of the plane
    def point_normal_to_plane(point, normal):
        a = normal.x
        b = normal.y
        c = normal.z
        d = -point.x * normal.x - point.y * normal.y - point.z * normal.z
        return Plane(a, b, c, d)

    # returns true if the plane contains the given point, else return false
    def contains(self, p):
        if (round(self.a * p.x + self.b * p.y + self.c * p.z + self.d,4) == 0):
            return True
        return False

    # returns a string the represents the plane
    def __str__(self):
        return f"{self.a} * x + {self.b} * y + {self.c} * z + {self.d} = 0"



class Segment:
  def __init__(self,p1_,p2_):
    self.p1 = p1_
    self.p2 = p2_




# global function Gets 2 points
# Calculates the plane that is perpendicular to the section connecting the two points and passing through the midpoint
def get_bisection_plane(p1, p2):
    mid = (p1 + p2).div(2)
    nvector = p2 - p1
    return Plane.point_normal_to_plane(mid,nvector)


# calculates an intersection point of 3 planes
def intersection_3Planes(pl1, pl2, pl3):
  A = np.array([[pl1.a, pl1.b, pl1.c], [pl2.a, pl2.b, pl2.c], [pl3.a, pl3.b, pl3.c]])
  b = np.array([-pl1.d, -pl2.d, -pl3.d])
  try:
    x = np.linalg.solve(A, b)
    p = Point(x[0], x[1], x[2])
  except:
    p = None
  finally:
    return p


# class that contains a point and all the middle vertical planes passing through the point
class VoronoiPoint:
   def __init__(self, point, planes):
     self.point = point
     self.planes = planes


# class that contains a segment and all the middle vertical planes passing through the segment
class VoronoiSegment:
   def __init__(self, segment, planes):
     self.segment = segment
     self.planes = planes


# a class that has a graph representation of a voronoi cell
class Graph:
    # gets the points, segemnts and bisections (between the main point and the points in the lattice)
    def __init__(self, points, segments, bisections):
        self.vertices = []
        self.edges = []
        self.faces = []

        vertice_num = dict()

        # representation of the vertices
        for i, p in enumerate(points):
            vertice_num["(" + str(round((p.point.x),4)) + "," + str(round((p.point.y),4)) + "," + str(round((p.point.z),4))+ ")"] = i
            self.vertices.append((str(i),len(p.planes)))

        # representation of the edges
        for i, s in enumerate(segments):
            p1 = "(" + str(round(s.segment.p1.x, 4)) + "," + str(round(s.segment.p1.y, 4)) + "," + str(round(s.segment.p1.z, 4))+ ")"
            p2 = "(" + str(round(s.segment.p2.x, 4)) + "," + str(round(s.segment.p2.y, 4)) + "," + str(round(s.segment.p2.z, 4))+ ")"
            str_s = "(" + str(vertice_num[p1]) + "," + str(vertice_num[p2]) + ")"
            self.edges.append((str_s,len(s.planes)))

        # representation of the faces
        for i,pl in enumerate(bisections):
            str_f = "("
            ps = []
            for p in points:
                if pl.contains(p.point):
                    ps.append(p)
            if len(ps) < 3:
                continue

            str_p = "(" + str(round(ps[0].point.x,4)) + "," + str(round(ps[0].point.y,4)) + "," + str(round(ps[0].point.z,4))+ ")"
            str_f += str(vertice_num[str_p])

            p = ps[0]
            check = 0
            segments_ = segments.copy()
            while check == 0 or p != (ps[0]).point:
                check += 1
                try:
                    p_ = p.point
                except:
                    p_ = p
                neighbor, segments_ = find_neighbor(segments_,p_,ps)
                if neighbor == ps[0].point:
                    break
                if neighbor != None:
                    str_neighbor = "(" + str(round(neighbor.x,4)) + "," + str(round(neighbor.y,4))+ "," + str(round(neighbor.z,4)) + ")"
                    str_f += "," + str(vertice_num[str_neighbor])
                    p = neighbor
                if neighbor == None:
                    break

            str_f += ")"
            self.faces.append(str_f)


# gets a list of segments, a point and we want to find its neighbor and a list of possible points
def find_neighbor(segments_,p,ps):
    neighbor = None
    for s in segments_:
        if s.segment.p1 == p:
            for vp in ps:
                if vp.point == s.segment.p2:
                    neighbor = s.segment.p2
                    segments_.remove(s)
                    return neighbor, segments_

        elif s.segment.p2 == p:
            for vp in ps:
                if vp.point == s.segment.p1:
                    neighbor = s.segment.p1
                    segments_.remove(s)
                    return neighbor, segments_
        else:
            continue
    return neighbor,segments_


# a class of the voronoi cell contains the points, segnents and the graph represntation of the cell.
class VoronoiCell:
    def __init__(self, points, segments, bisections):
        self.points = []
        self.segments = []
        for p in points:
            planes = []
            for pl in bisections:
                if pl.contains(p):
                    planes.append(pl)
            self.points.append(VoronoiPoint(p, planes))

        for s in segments:
            planes = []
            for pl in bisections:
                if pl.contains(s.p1) and pl.contains(s.p2):
                    planes.append(pl)
            self.segments.append(VoronoiSegment(s, planes))

        graph = Graph(self.points, (self.segments).copy(), bisections)
        self.graph_json = json.dumps(graph.__dict__)
        print(self.graph_json)


# creates the neighbors of the main point in the lattice by the spanning vectors
def neighbors_from_vectors(vector1, vector2, vector3, movings):
    points = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for l in [-1, 0, 1]:
                if i == 0 and j == 0 and l == 0:
                    continue
                p = Point.from_vector(vector1.mul(i) + vector2.mul(j) + vector3.mul(l))
                if movings == 'true':
                    p.x += float(random.uniform(-0.05, 0.05))
                    p.y += float(random.uniform(-0.05, 0.05))
                    p.z += float(random.uniform(-0.05, 0.05))
                points.append(p)

    ind = random.randint(0, len(points) -1)

    return points


# calculates all the intersection points between all triples in yhe given list of planes
def points_intersections_from_planes(planes):
    points_intersctions = []
    for i,pl1 in enumerate(planes):
      for j,pl2 in enumerate(planes[i+1:]):
        for pl3 in planes[i+j+2:]:
          p = intersection_3Planes(pl1,pl2,pl3)
          if p != None:
            points_intersctions.append(p)

    already_in = []
    points_intersctions_copy = points_intersctions.copy()
    points_intersctions = []
    for p in points_intersctions_copy:
        if (p.x, p.y, p.z) not in already_in:
            points_intersctions.append(p)
            already_in.append((p.x, p.y, p.z))

    return points_intersctions


# return 2 planes passing through the segement
def planes_from_segment(s):
    d = s.p2
    planes = []

    if round(d.x,4) == 0 and round(d.y,4) == 0:
        planes = [Plane(1,0,0,-d.x),Plane(0,1,0,-d.y)]
    elif round(d.y,4) == 0 and round(d.z,4) == 0:
        planes = [Plane(0, 0, 1, -d.z), Plane(0, 1, 0, -d.y)]
    elif round(d.x,4) == 0 and round(d.z,4) == 0:
        planes = [Plane(0, 0, 1, -d.z), Plane(1, 0, 0, -d.x)]
    elif round(d.x,4) == 0:
        planes = [Plane(1, 0, 0, -d.x), Plane(0,d.z,-d.y,0)]
    elif round(d.y,4) == 0:
        planes = [Plane(0, 1, 0, -d.y), Plane(d.z,0,-d.x,0)]
    elif round(d.z,4) == 0:
        planes = [Plane(0, 0, 1, -d.z), Plane(d.y,-d.x,0,0)]
    else:
      planes = [Plane(d.y,-d.x,0,0),Plane(0,d.z,-d.y,0)]
    return planes


# calculates which of the intersection points of the bisection are vertices of the voronoi cell
def get_voronoi_points(points_intersctions, bisections):
    segments = []
    for p in points_intersctions:
        segments.append(Segment(Point(0, 0, 0), p))

    voronoi_points = []
    not_voronoi_points = []

    i=0
    for s in segments:
        i += 1
        planes_of_s = planes_from_segment(s)
        j=0
        for pl in bisections:
            j += 1
            inter = intersection_3Planes(pl, planes_of_s[0], planes_of_s[1])

            if inter != None and not (((round(s.p2.x,4) >= 0 and (round(inter.x,4) >= round(s.p2.x,4) or round(inter.x,4) <= 0)) or (
            round(s.p2.x,4) <= 0 and (round(inter.x,4) <= round(s.p2.x,4) or round(inter.x,4) >= 0))) and (
            (round(s.p2.y,4) >= 0 and (round(inter.y,4) >= round(s.p2.y,4) or round(inter.y,4) <= 0)) or (
            round(s.p2.y,4) <= 0 and (round(inter.y,4) <= round(s.p2.y,4) or round(inter.y,4) >= 0))) and (
            (round(s.p2.z,4) >= 0 and (round(inter.z,4) >= round(s.p2.z,4) or round(inter.z,4) <= 0)) or (
            round(s.p2.z,4) <= 0 and (round(inter.z,4) <= round(s.p2.z,4) or round(inter.z,4) >= 0)))):

                not_voronoi_points.append(s.p2)
                break

    for p in points_intersctions:
        if p not in not_voronoi_points:
            voronoi_points.append(p)

    return voronoi_points


# gets the points of the cell and the possible planes of the cell, returns the segments of the cell
def get_segments(planes, points):
    segments = []
    for i, pl1 in enumerate(planes):
        for pl2 in planes[i + 1:]:
            pts = []
            for p in points:
                if pl1.contains(p) and pl2.contains(p):
                    pts.append(p)
            if len(pts) == 2:
                segments.append(Segment(pts[0], pts[1]))
    un = []
    for s in segments:
        for s2 in un:
            if (s.p1 == s2.p1 and s.p2 == s2.p2) or (s.p1 == s2.p2 and s.p2 == s2.p1):
                break
        else:
            un.append(s)
    return un


# get the voronoi cell by 3 spanning vectors of the lattice and a possible addition of movings
def get_cell(vector1, vector2, vector3, movings):

    points = neighbors_from_vectors(vector1, vector2, vector3, movings)

    bisections = []
    for p in points:
        bisections.append(get_bisection_plane(Point(0, 0, 0), p))

    points_intersctions = points_intersections_from_planes(bisections)

    voronoi_points = get_voronoi_points(points_intersctions, bisections)

    segments = get_segments(bisections, voronoi_points)

    cell = VoronoiCell(voronoi_points, segments, bisections)
    return cell

# get the voronoi cell by list of points
def get_cell_from_points(points):


    bisections = []
    for p in points:
        bisections.append(get_bisection_plane(Point(0, 0, 0), p))

    points_intersctions = points_intersections_from_planes(bisections)


    voronoi_points = get_voronoi_points(points_intersctions, bisections)

    segments = get_segments(bisections, voronoi_points)
    

    cell = VoronoiCell(voronoi_points, segments, bisections)
    return cell

# Visual representation in matplotlib
"""
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []
    cs = []
    colors = ['w', 'k', 'g', 'r', 'y', 'm', 'b', 'peru', 'lime', 'orange']
    for p in cell.points:
        xs.append(float(p.point.x))
        ys.append(float(p.point.y))
        zs.append(float(p.point.z))
        cs.append(colors[len(p.planes)])

    scatter = ax.scatter(xs, ys, zs, c =cs)
    ax.scatter([0], [0], [0])

    for l in cell.segments:
        ax.plot([l.segment.p1.x, l.segment.p2.x],[l.segment.p1.y, l.segment.p2.y], [l.segment.p1.z, l.segment.p2.z], colors[len(l.planes)])

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=str(i),
                            markerfacecolor=c, markersize=15)for i,c in enumerate(colors) ]
    ax.legend(handles=legend_elements, loc ='upper left')
    plt.show()"""


class VoronoiCellJson:
   def __init__(self, voronoicell):
    self.points = []
    for p in voronoicell.points:
        point = p.point
        self.points.append(((float(point.x),float(point.y),float(point.z)),len(p.planes)))

    self.segments = []
    for s in voronoicell.segments:
        point1 = s.segment.p1
        point2 = s.segment.p2
        self.segments.append((((float(point1.x),float(point1.y),float(point1.z)),(float(point2.x),float(point2.y),float(point2.z))),len(s.planes)))


def main(argv):
    if argv[1] == 'points':
        points = []
        for arg in argv[2:]:
            arr = arg.split(',')
            points.append(Point(float(arr[0]), float(arr[1]), float(arr[2])))
        cell = get_cell_from_points(points)
    else:
        vector1 = Vector(int(argv[1]),int(argv[2]),int(argv[3]))
        vector2 = Vector(int(argv[4]),int(argv[5]),int(argv[6]))
        vector3 = Vector(int(argv[7]),int(argv[8]),int(argv[9]))
        cell = get_cell(vector1, vector2, vector3, argv[10])
    cell_json = VoronoiCellJson(cell)

    jsonStr = json.dumps(cell_json.__dict__)
    print(jsonStr)

    return cell


if __name__ == '__main__':
    main(sys.argv)

