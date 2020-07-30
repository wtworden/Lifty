


from lifty.sage_types import *
import random
import itertools as it


COLORS = ['blue','red','green','cyan','black','grey','pink','purple','orange','brown','limegreen','olive']



class Triangulation:
    def __init__(self, edges, top_poly):
        for i in range(len(edges)):
            assert edges[i].index() == i
        self._edges = edges
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

        for i in range(len(self._vertices)):
            self._vertices[i]._index = i

        for i in range(len(self._vertices)):
            self.set_incident_edges(i)

        self._is_lifted = False

        self._linearized = None

    def set_incident_edges(self, vert_indx):
        v = self._vertices[vert_indx]
        p = v.point()
        L = [[e,e.point(1)] for e in self._edges if e.vert(0) == v] + [[e,e.point(-2)] for e in self._edges if e.vert(1) == v]
            ### TO DO: need to make sure the intervals give angles between adjacent edges that are intervals not containing zero.
            ###        otherwise, increase complexity. At this point we know all subvertices have disjoint intervals, and self is embedded.
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
        self._vertices[vert_indx]._incident_edges = incident_edges

    def top_poly(self):
        return self._top_poly

    def vertex(self,i):
        return self._vertices[i]

    def vertices(self):
        return self._vertices

    def num_vertices(self):
        return len(self._vertices)

    def edge(self,i):
        return self._edges[i]

    def edges(self):
        return self._edges

    def num_edges(self):
        return len(self._edges)

    def is_linear(self):
        for e in self.edges():
            if len(e.segments()) > 1:
                return False
        return True

    def linearized(self):
        if not self.is_linear():
            if self._linearized == None:
                straight_edges = []
                for i in range(self.num_edges()):
                    e = self.edge(i)
                    if self.is_lifted_tri():
                        straight_edges.append(MarkedEdge([e.point(0),e.point(-1)], e.maps_to(), e.index()))
                    else:
                        straight_edges.append(Edge([e.point(0),e.point(-1)], e.index()))
                if self.is_lifted_tri():
                    T = LiftedTriangulation(straight_edges, self.polynomial())
                else:
                    T = Triangulation(straight_edges)
                if not T.is_embedded():
                    return 'not embedded'
                for i in range(self.num_vertices()):
                    inds1 = self.vertex(i).incident_edges()
                    inds2 = T.vertex(i).incident_edges()
                    l = len(inds1)
                    if not inds1 in [[inds2[(i+j)%l] for i in range(l)] for j in range(l)]:
                        return 'vertex bad'
                self._linearized = T
            return self._linearized
        else:
            return None

    def add_edges_to_infinity(self):
        if not self.is_lifted_tri():
            for v in self.vertices():
                v_edges = v.incident_edges()
                l = len(v_edges)
                for i in range(l):
                    e1 = self.edge(v_edges[i])
                    e2 = self.edge(v_edges[(i+1)%l])
                    
                    

    def combinatorialized(self):
        pass

    def is_embedded(self):
        E = self._edges
        for i in range(self.num_edges()):
            e0 = self.edge(i)
            z0 = (e0.point(0).cif()+e0.point(-1).cif())/2
            r0 = abs(z0-e0.point(0).cif())
            for j in range(i+1,self.num_edges()):
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


    def plot(self, thickness=.7, DPI=400, show_vertices=True):
        G = Graphics()
        lines = []
        for i in range(self.num_edges()):
            edge = self.edge(i)
            points = edge.points()
            if self.is_lifted_tri():
                lines.append(Line([points[j].cif() for j in range(len(points))], color=COLORS[edge.maps_to()], thickness=thickness, zorder=10))
            else:
                lines.append(Line([points[j].cif() for j in range(len(points))], color=COLORS[i], thickness=thickness, zorder=10))
        for L in lines:
            G += L
        if show_vertices:
            for v in self._vertices:
                color = 'black'
                for c in self.top_poly().postcritical_set():
                    if v.point().cif().overlaps(c.cif()):
                        color = 'red'
                        break
                G += Pt(v.point().cif(), color=color, zorder=20)
        G.show(dpi=DPI)

class LiftedTriangulation(Triangulation):
    def __init__(self, marked_edges, top_poly):
        Triangulation.__init__(self, marked_edges, top_poly)

        self._is_lifted = True


class Edge:
    def __init__(self, points, index = None):
        self._vertices = [Vertex(points[0]),Vertex(points[-1])]
        self._points = points
        self._segments = [Segment([points[i],points[i+1]], self) for i in range(len(points)-1)]
        self._index = index
        self._triangulation = None

    def triangulation(self):
        return self._triangulation

    def vert(self,i):
        return self._vertices[i]

    def point(self,i):
        return self._points[i]

    def points(self):
        return self._points

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
    def __init__(self, points, maps_to, index):
        Edge.__init__(self, points, index)

        self._maps_to = maps_to

    def maps_to(self):
        return self._maps_to


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
    def __init__(self, point, index=None, incident_edges=None):

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

    def incident_edges(self):
        return self._incident_edges

    def __eq__(self,other):
        return self.point() == other.point()

    def __ne__(self,other):
        return self.point() != other.point()

class Point:
    def __init__(self, sage_alg_num, prec=53):
        self._alg = sage_alg_num
        self._cif = self._alg.interval(ComplexIntervalField(prec))
        self._precision = prec

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





