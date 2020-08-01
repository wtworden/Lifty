

from __future__ import print_function

from lifty.sage_types import *
import random
import itertools as it

from lifty.kernel.errors import IntervalError, IterationError, SubdivisionError


COLORS = ['blue','red','green','cyan','black','grey','pink','purple','orange','brown','limegreen','olive']



class Triangulation:
    def __init__(self, edges, top_poly):
        for i in range(len(edges)):
            assert edges[i].index() == i
        self._edges = edges
        self._num_finite_edges = len(edges)
        self._top_poly = top_poly

        for edge in edges:
            edge._triangulation = self

        self._vertices = []
        for e in self._edges:
            for i in range(2):
                v = e.vert(i)
                if v not in self._vertices:
                    v._triangulation = self
                    self._vertices.append(v)
        self._vertices.append(Vertex(Point(Infinity)))

        for i in range(len(self._vertices)):
            self._vertices[i]._index = i

        for i in range(len(self._vertices)-1):
            self.set_incident_edges(i)

        self._is_lifted = False

        self._linearized = None

    def set_incident_edges(self, vert_indx):
        v = self.vertices()[vert_indx]
        p = v.point()
        L = [[e,e.point(1)] for e in self.edges() if e.vert(0) == v] + [[e,e.point(-2)] for e in self.edges() if e.vert(1) == v]


        if len(L) > 1:
            # make sure the ordering of segments coming from v can be determined using intervals.
            for j in range(len(L)):
                q1, q2 = L[j][1], L[(j+1)%(len(L))][1]
                prec = max([q1.precision(), q2.precision()])
                arg1 = (q1.cif() - p.cif()).arg()
                arg2 = (q2.cif() - p.cif()).arg()
                while arg1.overlaps(arg2):
                    prec = 2*prec
                    if prec < self.top_poly().max_precision():
                        q1.set_precision(prec)
                        q2.set_precision(prec)
                    else:
                        raise IntervalError('max precision reached: cyclic ordering of edges around vertex {} cannot be determined.'.format(vert_indx))

        L.sort(key = lambda x: (x[1].cif()-p.cif()).arg())
        incident_edges = [e[0].index() for e in L]
        self.vertices()[vert_indx]._incident_edges = incident_edges

    def top_poly(self):
        return self._top_poly

    def vertex(self,i):
        return self._vertices[i]

    def vertices(self):
        return self._vertices

    def finite_vertices(self):
        return self._vertices[:-1]

    def num_vertices(self):
        return len(self._vertices)

    def num_finite_vertices(self):
        return len(self._vertices) - 1

    def edge(self,i):
        return self._edges[i]

    def edges(self):
        return self._edges

    def finite_edges(self):
        return self._edges[:self.num_finite_edges()]

    def num_finite_edges(self):
        return self._num_finite_edges

    def num_infinite_edges(self):
        return self.num_edges() - self.num_finite_edges()

    def num_edges(self):
        return len(self._edges)

    def is_linear(self):
        for e in self.finite_edges():
            if len(e.segments()) > 1:
                return False
        return True

    def add_edges_to_infinity(self):
        if self.num_infinite_edges() == 0: # only do something if we don't already have edges to infty
            if not self.is_lifted_tri():
                verts = self.finite_vertices()
                inf_vert = self.num_vertices() - 1

                # for each vertex we determine if its on the boundary (since the pc_tri is convex, we just look
                # for an angle between incident edges that is >= pi), then if it is we put the edge to infinity halfway 
                # between (in angle) these two edges.
                for v in verts:
                    L = len(v.incident_edges())
                    edges_to_infty = []

                    # for each i, check if the angle between incident edge i and i+1 (mod L) is >= pi
                    for i in range(L):
                        e1, e2 = v.incident_edges()[i], v.incident_edges()[(i+1)%L]
                        v1, v2 = self.edge(e1).other_vert(v), self.edge(e2).other_vert(v)
                        angle = ((v2.point().cif() - v.point().cif())*exp(-I*(v1.point().cif() - v.point().cif()).arg())).arg()
                        print(angle)
                        # sage returns arg in the range [-pi,pi], so an angle between e1 and e2 of >= pi
                        # corresponds here to angle <= 0 or angle = pi
                        a1 = v1.point().alg()
                        a2 = v2.point().alg()
                        a = v.point().alg()
                        print((a1,a2,a))
                        if angle.overlaps(RIF(pi)):
                            b = (a1-a).imag() - QQbar(I)*(a1-a).real()
                        elif angle.overlaps(RIF(0)):
                            b = -(a1-a)
                        elif angle < 0:
                            b = -((a2-a)/(a2-a).abs() + (a1-a)/(a1-a).abs())/2

                        s = Point((b/b.abs()) + a)
                        e = Edge([v,self.vertices()[-1]],[v.point(), s, Point(Infinity)])
                        self._edges.append(e)
                        e._index = self.num_edges()-1
                        edges_to_infty.append((e,(i+1)%L))
                    edges_to_infty.sort(key=lambda x: x[1])
                    edges_to_infty.reverse()
                    for e,j in edges_to_infty:
                        v.incident_edges().insert(j,e.index())
                E = self.edge(self.num_finite_edges())
                inf_edges_ordered = [E.index()]
                v = E.vert(0)
                while len(inf_edges_ordered) < self.num_infinite_edges():
                    next_index = (v.incident_edges().index(E.index()) + 1)%(v.valence())
                    e = self.edge(v.incident_edges()[next_index])
                    v = e.other_vert(v)
                    next_index = (v.incident_edges().index(e.index()) + 1)%(v.valence())
                    E = self.edge(v.incident_edges()[next_index])
                    inf_edges_ordered.append(E.index())
                inf_edges_ordered.reverse()
                self.vertex(-1)._incident_edges = inf_edges_ordered
            else:
                pc_tri = self.top_poly().pc_triangulation()
                pc_tri.add_edges_to_infinity()



                    

    def combinatorialized(self):
        pass

    def is_embedded(self):
        E = self._edges
        for i in range(self.num_finite_edges()):
            e0 = self.edge(i)
            z0 = (e0.point(0).cif()+e0.point(-1).cif())/2
            r0 = abs(z0-e0.point(0).cif())
            for j in range(i+1,self.num_finite_edges()):
                e1 = self.edge(j)
                z1 = (e1.point(0).cif()+e1.point(-1).cif())/2
                r1 = abs(z1-e1.point(0).cif())
                if abs(z0 - z1) > r0 + r1:
                    continue
                else:
                    for k in range(e1.num_segments()):
                        seg1 = e1.segment(k)
                        x1 = (seg1.endpoint(0).cif() + seg1.endpoint(1).cif())/2
                        s1 = abs(x1 - seg1.endpoint(0).cif())
                        if abs(z0 - x1) > r0 + s1:
                            continue
                        else:
                            for m in range(e0.num_segments()):
                                seg0 = e0.segment(m)
                                x0 = (seg0.endpoint(0).cif() + seg0.endpoint(1).cif())/2
                                s0 = abs(x0 - seg0.endpoint(0).cif())
                                if abs(x0 - x1) > s0 + s1:
                                    continue    
                                else:
                                    if seg0.transverse_to(seg1):
                                        return False, [(i,m),(j,k)]
                                    else:
                                        continue
        return True, None

    def subdivide_segment(self, edge_indx, seg_indx, subdivision_point):
        e = self._edges[edge_indx]
        e.subdivide_segment(seg_indx, subdivision_point)

        if seg_indx == 0:
            v = self.edge(edge_indx).vert(0)
            self.set_incident_edges(v.index())

        elif seg_indx == self.edge(edge_indx).num_segments() - 2:
            v = self.edge(edge_indx).vert(1)
            self.set_incident_edges(v.index())


    def is_lifted_tri(self):
        return self._is_lifted


    def plot(self, thickness=.7, DPI=400, aspect_ratio='automatic',show_vertices=True):
        G = Graphics()
        lines = []
        for i in range(self.num_edges()):
            edge = self.edge(i)
            points = edge.finite_points()
            ### for some reason Sage will not plot a line if its endpoints are elements of the ComplexIntervalField
            ### with imaginary part =0. So, as a somewhat hacky workaround we'll just add a very small multiple of I
            ### to each point.
            eps = CIF(0.00000000001*I)
            if self.is_lifted_tri():
                lines.append(Line([points[j].cif()+eps for j in range(len(points))], color=COLORS[edge.maps_to()], thickness=thickness, zorder=20))
            else:
                f = self.num_finite_edges()
                if i < f:
                    color = COLORS[i%12]
                else:
                    color = COLORS[f]
                lines.append(Line([points[j].cif()+eps for j in range(len(points))], color=color, thickness=thickness, zorder=20))
        for L in lines:
            G += L
        if show_vertices:
            for v in self.finite_vertices():
                color = 'black'
                for c in self.top_poly().postcritical_set():
                    if v.point().cif().overlaps(c.cif()):
                        color = 'red'
                        break
                G += Pt(v.point().cif(), color=color, zorder=30)
        G.show(dpi=DPI, aspect_ratio=aspect_ratio)

class LiftedTriangulation(Triangulation):
    def __init__(self, marked_edges, top_poly):
        Triangulation.__init__(self, marked_edges, top_poly)

        self._pc_triangulation = top_poly.pc_triangulation()
        self._is_lifted = True
        self._vertices[-1] = MarkedVertex(Point(Infinity),self._pc_triangulation.vertices()[-1])

class Edge:
    def __init__(self, vertices, points, index = None):
        self._vertices = vertices
        self._points = points
        if points[-1].alg() == Infinity:
            self._finite_points = points[:-1]
            self._ray = Ray([points[-2], points[-1]], self)
        else:
            self._finite_points = points
            self._ray = None
        self._segments = [Segment([self._finite_points[i],self._finite_points[i+1]], self) for i in range(len(self._finite_points)-1)]
        self._index = index
        self._triangulation = None

    def triangulation(self):
        return self._triangulation

    def vert(self,i):
        return self._vertices[i]

    def other_vert(self, vert):
        if vert == self.vert(0):
            return self.vert(1)
        elif vert == self.vert(1):
            return self.vert(0)

    def point(self,i):
        return self._points[i]

    def points(self):
        return self._points

    def finite_points(self):
        return self._finite_points

    def segment(self,i):
        return self._segments[i]

    def num_segments(self):
        return len(self._segments)

    def segments(self):
        return self._segments

    def index(self):
        return self._index

    def subdivide_segment(self, seg_indx, subdivision_point):
        self._segments[seg_indx] = Segment([self.segment(seg_indx).endpoint(0),subdivision_point], self)
        self._segments.insert(seg_indx+1, Segment([subdivision_point, self.segment(seg_indx).endpoint(1)], self))

        self._points.insert(seg_indx+1, subdivision_point)


class MarkedEdge(Edge):
    def __init__(self, vertices, points, maps_to, index):
        Edge.__init__(self, vertices, points, index)

        self._maps_to = maps_to

    def maps_to(self):
        return self._maps_to

class Ray:
    def __init__(self, endpoints, edge):
        self._endpoints = endpoints
        self._edge = edge

    def edge(self):
        return self._edge

    def endpoint(self,i):
        return self._endpoints[i]

class Segment:
    def __init__(self, endpoints, edge):
        assert endpoints[0] != endpoints[1]
        self._endpoints = endpoints
        self._edge = edge
        self._slope = None

    def __eq__(self, other):
        if set(self._endpoints) == set(other._endpoints):
            return True
        else:
            return False

    def __ne__(self,other):
        if not self.__eq__(other):
            return True
        else:
            return False

    def edge(self):
        return self._edge

    def endpoint(self,i):
        return self._endpoints[i]

    def slope(self): ## return an RIF, not a algebraic integer
        a,b = self.endpoint(0).cif(), self.endpoint(1).cif()
        if self._slope == None:
            if a.real().overlaps(b.real()):
                self._slope = Infinity
            else:
                self._slope = (a.imag() - b.imag())/(a.real() - b.real())
        return self._slope

    def y_intercept(self): ## return an RIF, not a algebraic integer
        if self.slope() != Infinity:
            return -self.slope()*self.endpoint(0).cif().real() + self.endpoint(0).cif().imag()
        else:
            return None

    def transverse_to(self, other):
        a, b = self.endpoint(0), self.endpoint(1)
        c, d = other._endpoints[0], other._endpoints[1]
        if ( self.above_seg(c) and self.below_seg(d) ) or (self.above_seg(d) and self.below_seg(c) ):
            if ( other.above_seg(a) and other.below_seg(b) ) or (other.above_seg(b) and other.below_seg(a) ):
                return True
        return False

    def is_endpoint(self,point):
        a,b = self.endpoint(0).cif(), self.endpoint(1).cif()
        c = point.cif()
        if a.overlaps(c) or b.overlaps(c):
            return True
        else:
            return False



    def above_seg(self, point):
        a,b,c = self.endpoint(0).cif(), self.endpoint(1).cif(), point.cif()
        s = self.slope()
        if self.slope() == Infinity:
            if c.real() > a.real().union(b.real()):
                return True
        elif self.slope().overlaps(RIF(0)):
            if c.imag() > a.imag().union(b.imag()):
                return True
        else:
            if c.imag() > self.slope()*c.real() + self.y_intercept():
                return True
        return False

    def below_seg(self, point):
        a,b,c = self.endpoint(0).cif(), self.endpoint(1).cif(), point.cif()
        s = self.slope()
        if self.slope() == Infinity:
            if c.real() < a.real().union(b.real()):
                return True
        elif self.slope().overlaps(RIF(0)):
            if c.imag() < a.imag().union(b.imag()):
                return True
        else:
            if c.imag() < self.slope()*c.real() + self.y_intercept():
                return True
        return False

    def on_seg(self, point):
        a,b,c = self.endpoint(0).cif(), self.endpoint(1).cif(), point.cif()
        s = self.slope()
        if self.slope() == Infinity:
            if c.real().overlaps(a.real().union(b.real())):
                if c.imag().overlaps(a.imag().union(b.imag())):
                    return True
        elif self.slope().overlaps(RIF(0)):
            if c.imag().overlaps(a.imag().union(b.imag())):
                if c.real().overlaps(a.real().union(b.real())):
                    return True
        else:
            if c.imag().overlaps(self.slope()*c.real() + self.y_intercept()):
                if c.real().overlaps(a.real().union(b.real())):
                    return True
        return False


class Vertex:
    def __init__(self, point, index=None, incident_edges=[]):

        self._point = point
        self._index = index
        self._incident_edges = incident_edges
        self._triangulation = None

    def triangulation(self):
        return self._triangulation

    def point(self):
        return self._point

    def index(self):
        return self._index

    def valence(self):
        return len(self.incident_edges())

    def incident_edges(self):
        return self._incident_edges

    def __eq__(self,other):
        return self.point() == other.point()

    def __ne__(self,other):
        return self.point() != other.point()

class MarkedVertex(Vertex):
    def __init__(self, point, maps_to, index=None, incident_edges=[]):
        Vertex.__init__(self, point, index, incident_edges)

        self._maps_to = maps_to

    def maps_to(self):
        return self._maps_to


class Point:
    def __init__(self, sage_alg_num, prec=53):
        self._alg = sage_alg_num
        if self._alg == Infinity:
            self._cif = CIF(Infinity)
            self._is_infinity = True
        else:
            self._cif = self._alg.interval(ComplexIntervalField(prec))
            self._is_infinity = False
        self._precision = prec

    def is_infinity(self):
        return self._is_infinity

    def alg(self):
        return self._alg

    def cif(self):
        return self._cif

    def precision(self):
        return self._precision

    def set_precision(self, prec):
        self._precision = prec
        self._cif = self._alg.interval(ComplexIntervalField(prec))

    def overlaps(self,other):
        return self.cif().overlaps(other.cif())

    def __eq__(self,other):
        return type(self) == type(other) and self.cif().overlaps(other.cif())

    def __ne__(self,other):
        return not self.__eq__(other) 

    def __repr__(self):
        return self.cif().__repr__()

    def __str__(self):
        return self.cif().__str__()

    def __mul__(self,other):
        return Point(self.alg()*other)

    def __add__(self,other):
        return Point(self.alg() - other)

    def __add__(self,other):
        return Point(self.alg() - other)





