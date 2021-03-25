

from __future__ import print_function

from lifty.sage_types import *
import random
import itertools as it
from copy import deepcopy
#from curver import create_triangulation
#from curver.kernel.triangulation import Triangulation as curverTriangulation
#from curver.kernel.triangulation import Edge as curverEdge
#from curver.kernel.triangulation import Triangle as curverTriangle

from lifty.kernel.combinatorial_tri import *

from lifty.kernel.errors import IntervalError, IterationError, SubdivisionError


COLORS = ['blue','red','green','cyan','black','pink','purple','orange','brown','limegreen','olive','yellow','lightgrey']



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

        len_pc = len(self._top_poly.postcritical_set())
        i=len_pc
        for v in self._vertices:
            if v.point() in self._top_poly.postcritical_set():
                v._index = self._top_poly.postcritical_set().index(v.point())
            else:
                v._index = i
                i += 1
        self._vertices.sort(key = lambda v: v._index)
        self._vertices.append(Vertex(Point(Infinity)))

        for i in range(len(self._vertices)):
            self._vertices[i]._index = i

        for i in range(len(self._vertices)-1):
            self.set_incident_edges(i)

        self._is_lifted = False
        self._combinatorial_old = None
        self._combinatorial = None

    def set_incident_edges(self, vert_indx):
        v = self.vertices()[vert_indx]
        p = v.point()
        L = [[e,e.point(1)] for e in self.edges() if e.vert(0) == v] + [[e,e.point(-2)] for e in self.edges() if e.vert(1) == v]


        if len(L) > 1:
            # make sure the ordering of segments coming from v can be determined using intervals.
            for i in range(len(L)-1):
                for j in range(i+1,len(L)):
                    q1, q2 = L[i][1], L[j][1]
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
        # sort the edges so they are order counter-clockwise around v.
        L.sort(key = lambda x: (x[1].cif()-p.cif()).arg())
        incident_edges = [e[0].index() for e in L]
        self.vertices()[vert_indx]._incident_edges = incident_edges

    def min_vertex_angle(self):
        min_angle = 2*pi
        for v in self.vertices()[:-1]:
            p = v.point()
            incident = v.incident_edges()
            for i in range(len(incident)):
                e1, e2 = self.edge(incident[i]), self.edge(incident[(i+1)%len(incident)])
                p1 = e1.point(1) if p == e1.point(0) else e1.point(-2)
                p2 = e2.point(1) if p == e2.point(0) else e2.point(-2)
                angle = abs((p1.cif() - p.cif()).arg() - (p2.cif() - p.cif()).arg())
                min_angle = min(min_angle, angle)
        return min_angle

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

    def infinite_edges(self):
        return self._edges[self.num_finite_edges():]

    def num_edges(self):
        return len(self._edges)

    def is_linear(self):
        for e in self.finite_edges():
            if len(e.segments()) > 1:
                return False
        return True

    def add_edges_to_infinity(self):
        if self.num_infinite_edges() == 0: # only do something if we don't already have edges to infty
            verts = self.finite_vertices()
            inf_vert = self.vertex(-1)
            if not self.is_lifted_tri():

                # for each vertex v we determine if its on the boundary (since the pc_tri is convex, we just look
                # for an angle between incident edges that is >= pi), then if it is we put an edge to infinity
                # from that point so that the edge is parallel to one of the edges coming in to that point (unless
                # the angle is pi, in which case it may be orthogonal to the two incoming edges). We construct the edge 
                # this way because it contains the point (denoted s below) 2*a-a1, where a is the point at v,
                # and a1 is the point at one of the vertices that shares an edge with v. The point s constructed this
                # way appears to be much easier for Sage to lift via P (the topological polynomial) than some other 
                # reasonable choices.

                def on_perimeter(w):
                    L = len(w.incident_edges())
                    for i in range(L):
                        E1, E2 = w.incident_edges()[i], w.incident_edges()[(i+1)%L]
                        v1, v2 = self.edge(E1).other_vert(w), self.edge(E2).other_vert(w)
                        angle = ((v2.point().cif() - w.point().cif())*exp(-I*(v1.point().cif() - w.point().cif()).arg())).arg()
                        if angle.overlaps(RIF(0)) or angle < 0:
                            return True, (w, angle, E1, E2)

                # first we find a vertex on the perimeter of the finite part of the triangulation that has 
                # angle > pi between two consecutive incident edges (note the strict inequality---we don't
                # want the vertex we start with to be one with angle = pi incident edges).

                for w in verts:
                    boo, tup = on_perimeter(w)
                    if boo:
                        break

                lifted_pcs_center = sum([p.cif() for p in self.top_poly().lifted_pcs()])/len(self.top_poly().lifted_pcs())
                lifted_pcs_radius = max([abs(lifted_pcs_center-p.cif()) for p in self.top_poly().lifted_pcs()])


    #                                                   o  s = 2a - a1
    #                                                    \
    #                                                     \
    #                                                      \ e
    #                                                       \
    #                                                        \   ___ angle
    #                                                          /     \
    #                            a2              E2           | \ a   |
    #                               o ________________________V_ o    |
    #                                                             \ _/
    #                                                              \
    #                                                               \
    #                                                                \  E1
    #                                                                 \
    #                                                                  \
    #                                                                   o a1
    #                                                                      
    #                                ^
                ## see diagram above |

                w, angle, E1, E2 = tup
                v1 = self.edge(E1).other_vert(w)
                v2 = self.edge(E2).other_vert(w)
                a1 = v1.point().alg()
                a = w.point().alg()
                b = a - a1
                s = Point(a + b)

                # we want s and all of its lifts via P to be outside of a circle containing all 
                # lifts of the postcritical set. This will insure that edges to infinity of both the 
                # post-critical triangulation and the lifted triangulation are not truncated until
                # they are outside of this circle, and hence we do not lose any important intersection
                # information.

                def extend_edge(s, a, a1, b):

                    s_lifts = self.top_poly().lifts(s)
                    s_lifts.append(s)
                    inside = [abs(l.cif() - lifted_pcs_center) <= lifted_pcs_radius for l in s_lifts]
                    while any(inside):
                        b += a - a1
                        s = Point(a + b)
                        s_lifts = self.top_poly().lifts(s)
                        s_lifts.append(s)
                        inside = [abs(l.cif() - lifted_pcs_center) <= lifted_pcs_radius for l in s_lifts]
                    return s

                s = extend_edge(s, a, a1, b)

                e = Edge([w, inf_vert],[w.point(), s, inf_vert.point()])
                self._edges.append(e)
                e._index = self.num_edges()-1
                L = len(w.incident_edges())
                E1_index = w.incident_edges().index(E1)
                w.incident_edges().insert((E1_index+1)%L,e.index())

                prev_angle = angle
                prev_b = b

                v = v2
                e1 = E2
                e1_index = v.incident_edges().index(e1)
                L = len(v.incident_edges())
                e2 = v.incident_edges()[(e1_index+1)%L]

                # walk counter-clockwise around the perimeter of the finite triangulation, adding infinite 
                # edges as we go. Some vertices might appear on the perimeter more than one (e.g., for the 
                # airplane), so we need to keep going until we see the vertex w with previous and next edges
                # E1 and E2, resp.
                while (v,e1) != (w,E1):
                    L = len(v.incident_edges())
                    edges_to_infty = []
 

                    v1, v2 = self.edge(e1).other_vert(v), self.edge(e2).other_vert(v)
                    angle = ((v2.point().cif() - v.point().cif())*exp(-I*(v1.point().cif() - v.point().cif()).arg())).arg()

                    # sage returns arg in the range [-pi,pi], so an angle between e1 and e2 of >= pi
                    # corresponds here to angle <= 0 or angle = pi
                    assert angle.overlaps(RIF(pi)) or angle.overlaps(RIF(0)) or angle < 0
                    a1 = v1.point().alg()
                    a = v.point().alg()
                    if angle.overlaps(RIF(pi)):
                        if prev_angle.overlaps(RIF(0)):
                            b = (a1 - a)*QQbar(I)
                        else:
                            b = prev_b
                    else:
                        b = a - a1
                    s = Point(a + b)

                    s = extend_edge(s, a, a1, b)

                    e = Edge([v, inf_vert],[v.point(), s, inf_vert.point()])
                    self._edges.append(e)
                    e._index = self.num_edges()-1
                    v.incident_edges().insert((e1_index+1)%L,e.index())

                    prev_angle = angle
                    prev_b = b
                    v = v2
                    e1 = e2
                    e1_index = v.incident_edges().index(e1)
                    L = len(v.incident_edges())
                    e2 = v.incident_edges()[(e1_index+1)%L]


            else:
                pc_tri = self.top_poly().pc_triangulation()
                pc_tri.add_edges_to_infinity()
                P = self.top_poly()
                next_edge_index = self.num_edges()

                vertex_pts = [v.point() for v in self.vertices()]
                for e in pc_tri.infinite_edges():
                    lifts = P.lift_edge(e)
                    print('edge {} lifted.'.format(e.index()))
                    for lift in lifts:
                        j = vertex_pts.index(lift[0])
                        v,w = self.vertices()[j], self.vertices()[-1]
                        lift[0] = v.point()
                        lift[-1] = w.point()
                        edge = MarkedEdge([v,w],lift,e.index(),next_edge_index)
                        self._edges.append(edge)
                        next_edge_index += 1
                for i in range(len(self.finite_vertices())):
                    self.set_incident_edges(i)
#               F = pc_tri.num_finite_edges()
#               for v in verts:
#                   m = len(v.incident_edges())
#                   insertions = []
#                   for i in range(m):
#                       ind1, ind2 = v.incident_edges()[i], v.incident_edges()[(i+1)%m]
#                       e1, e2 = self.edge(ind1), self.edge(ind2)
#                       V = v.maps_to()
#                       IND1, IND2 = e1.maps_to(), e2.maps_to()
#                       j = V.incident_edges().index(IND1)
#                       M = len(V.incident_edges())
#                       next_edge_index = V.incident_edges()[(j+1)%M]
#                       if next_edge_index >= F:
#                           insertions.append(((i+1)%m, next_edge_index))
#                       else:
#                           assert next_edge_index == IND2
#                   insertions.sort()
#                   insertions.reverse()
#                   for i, edge_ind in insertions:
#                       new_edge = MarkedEdge([v, inf_vert],[v.point(), inf_vert.point()], edge_ind)
#                       self._edges.append(new_edge)
#                       new_edge._index = self.num_edges()-1
#                       v.incident_edges().insert(i,new_edge.index())

            E = self.edge(self.num_finite_edges())
            inf_edges_ordered = [E.index()]
            v = E.vert(0)
            while len(inf_edges_ordered) < self.num_infinite_edges():
                next_index = (v.incident_edges().index(E.index()) + 1)%(v.degree())
                e = self.edge(v.incident_edges()[next_index])
                v = e.other_vert(v)
                next_index = (v.incident_edges().index(e.index()) + 1)%(v.degree())
                E = self.edge(v.incident_edges()[next_index])
                inf_edges_ordered.append(E.index())
            inf_edges_ordered.reverse()
            self.vertex(-1)._incident_edges = inf_edges_ordered

    def ribbon_graph(self):
        L = self.num_edges()
        rho = PermutationGroupElement([(i+1,i+1+L) for i in range(L)])
        sigma_list = []
        for v in self.vertices():
            incident = v.incident_edges()
            incident_darts = []
            for i in incident:
                e = self.edge(i)
                if e.vertex(0) == v:
                    incident_darts.append(i+1)
                elif e.vertex(1) == v:
                    incident_darts.append(i+1+L)
            sigma_list.append(tuple(incident_darts))
        sigma = PermutationGroupElement(sigma_list)
        return RibbonGraph(sigma,rho)


    def combinatorial(self):
        if self._combinatorial is None:
            triangles = []
            edges = {}
            c_vertices = [CVertex(i) for i in range(len(self.vertices()))]

            # start with an edge, then build the triangles adjacent to the positive and negative reps
            # of the edge, but only if the triangle hasn't already been built (i.e., if the label does
            # not already appear as a key in the dictionary edges)
            for e in self.edges():
                for label in [int(e.index()), ~int(e.index())]:
                    if label not in edges:
                        e0 = e
                        ind = e.index()

                        # record the sign of the label, and the e-index (next_vert) of vertex at the forward end
                        # of e w.r.t. counter-clock-wise travel around the perimeter of the triangle
                        if label >= 0:
                            next_vert = 1
                            tri_signs = [1]
                        else:
                            next_vert = 0
                            tri_signs = [-1]                        

                        # if the opposite sign label (~label) is already in edges, then we don't want to create a new edge,
                        # but rather just point label to the edge that ~label points to.
                        if ~label not in edges:
                            edges[label] = CEdge(ind,[c_vertices[e0.vertex(0).index()], c_vertices[e0.vertex(1).index()]])
                        else:
                            edges[label] = edges[~label]

                        # initiate a list of edges for the triangle with the edge e
                        tri_edges = [edges[label]]

                        for i in range(2):

                            # get the Vertex associated to the e-index next_vert
                            v0 = e0.vertex(next_vert)
                            d = v0.degree()

                            # get the edge incident on v0 that comes before e0 in the clockwise ordering of 
                            # edges around v0. This will be the next edge of the triangle after e0 as we traverse
                            # the boundary of the triangle counter-clock-wise.
                            prev = (v0.incident_edges().index(ind)-1)%d

                            # get the index of the edge, then the actual Edge object
                            ind = v0.incident_edges()[prev]
                            e0 = self.edge(ind)

                            # get the sign of this edge, and the next vertex
                            if e0.vertex(0) == v0:
                                label = ind
                                next_vert = 1
                                tri_signs.append(1)
                            elif e0.vertex(1) == v0:
                                label = ~int(ind)
                                next_vert = 0
                                tri_signs.append(-1)

                            # create the CEdge if we haven't already created it with ~label
                            if ~label not in edges:
                                edges[label] = CEdge(ind,[c_vertices[e0.vertex(0).index()], c_vertices[e0.vertex(1).index()]])
                            else:
                                edges[label] = edges[~label]
                            tri_edges.append(edges[label])

                        #get the next verted, then get the incident edge that will should correspond to the
                        # edge we started with
                        v0 = e0.vertex(next_vert)
                        d = v0.degree()
                        prev = (v0.incident_edges().index(ind)-1)%d
                        ind = v0.incident_edges()[prev]
                        assert ind == e.index()
                        triangle = CTriangle(tri_edges, tri_signs)
                        triangles.append(triangle)
            vertex_markings = {}
            for i in range(len(self.vertices())):
                v = self.vertices()[i]
                if v.point() in self.top_poly().postcritical_set():
                    vertex_markings[i] = [True, True]
                else:
                    vertex_markings[i] = [False, False]
                if  v.point().is_infinity():
                    vertex_markings[i][1] = True
            if self.is_lifted_tri():
                vertex_mappings = [v.maps_to().index() for v in self.vertices()]
            else:
                edge_mappings = None
                vertex_mappings = None
            self._combinatorial = CTriangulation(triangles, c_vertices, vertex_markings, vertex_mappings, self.is_lifted_tri())
        
        # encode the edges of the postcritical triangulation as a multi-arc in this combinatorial triangulation
        if self._is_lifted:
            pc_tri = self._pc_triangulation
            arcs = []
            self_copy = deepcopy(self)
            eps = QQbar(QQ(1)/(10**(int(-log(min([seg.length() for edge in self_copy.edges() for seg in edge.segments()]).center()/1000,10)))))
            min_angle = self_copy.min_vertex_angle()
            vertices = [[None,None] for i in range(pc_tri.num_edges())]
            for pc_edge in pc_tri.edges():
                assert pc_edge.num_segments() == 1
                for i in range(2):
                    a = pc_edge.vertex(i).point().alg()
                    for j in range(self.num_vertices()):
                        v = self.vertex(j)
                        if a == v.point().alg():
                            vertices[pc_edge.index()][i] = j
                pc_edge_seg = pc_edge.segment(0)

                # first we'll make self_copy transverse to pc_edge, in such a way that every segment of self_copy
                # is transverse to pc_edge (i.e., no endpoint of a segment lies on pc_edge, except of course at
                # the vertices of pc_edge).


                for e in self_copy.edges():
                    for i in range(len(e.points())):
                        p = e.point(i)
                        # if p in on pc_edge but is not at one of its vertices...
                        if pc_edge_seg.on_seg(p) and not ( p==pc_edge.vertex(0).point() or p==pc_edge.vertex(1).point() ):

                            # push p off of pc_edge_seg to the left (relative to orientation of pc_edge_seg) by eps. By
                            # definition eps is about .001 times the length of the sortest segment of self_copy,
                            # so the amount we push off is small compared to the length of any segment, but large
                            # compared to the precisions of intervals (which is less than .000001 times the length
                            # of the shortest segment). Before we do the pushoff, we need to make sure eps is small
                            # enough so that we stay away from the other edges of pc_tri. This can be done by
                            # considering the smallest angle min_angle of pc_tri, and ensuring that eps < x*tan(min_angle)/2,
                            # where x in the mininal distance from p to a vertex of pc_edge.

                            x = min([abs(p.cif() - pc_edge.vertex(0).point().cif()), abs(p.cif() - pc_edge.vertex(1).point().cif())])
                            y = x*tan(min_angle)/2
                            eps0 = eps
                            while eps0 > y:
                                eps0 = eps0/10
                                prec = min([point.precision() for edge in self_copy.edges() for point in edge.points()])
                                if eps0 < 1000*(10**(-prec*log(2,10))):
                                    if prec < self.top_poly().max_precision():
                                        prec = prec*2
                                        for point in [point for edge in self_copy.edges() for point in edge.points()]:
                                            point.set_precision(prec)
                                    else:
                                        raise IntervalError('max precision reached: when trying to compute combinatorial triangulation')
                            p._set_alg(p.alg()+eps0*pc_edge_seg.unit_tangent()*QQbar(I))

            alg_intersections_lists = []
            alg_arcs = [[0 for i in range(self_copy.num_edges())] for j in range(pc_tri.num_edges())]
            geom_arcs = [[0 for i in range(self_copy.num_edges())] for j in range(pc_tri.num_edges())]
            for lifted_edge in self_copy.edges():
                alg_intersections = []
                for i in range(len(lifted_edge.segments())):
                    seg = lifted_edge.segment(i)
                    for e in pc_tri.edges():
                        if seg.transverse_to(e.segment(0)):
                            a,b = seg.endpoint(0).cif(), seg.endpoint(1).cif()
                            c,d = e.segment(0).endpoint(0).cif(), e.segment(0).endpoint(1).cif()
                            # first translate by c to get c==0
                            a,b = a-c, b-c
                            c,d = c-c, d-c
                            # now rotate so that (c)-->--(d) has slope 1
                            z = exp(I*(pi/4 - d.arg()))
                            a,b = a*z, b*z
                            d = d*z
                            # seg is now parametrized by P+tv, with P=a and v=b-a.
                            v = b-a
                            P = a
                            # the parameter t at the intersection point is:
                            t = (P.imag()-P.real())/(v.real() - v.imag())
                            alg_int = 2*int((d.real()*(b-a).imag() - (b-a).real()*d.imag()) > 0) - 1
                            alg_intersections.append((i+t,alg_int,e.index()))
                alg_intersections.sort(key=lambda x:x[0])
                # now that the intersections are in order, we don't need the t+i parameter
                print(lifted_edge.index(),alg_intersections)
                alg_ints_tuples = [(tup[1],tup[2]) for tup in alg_intersections]
                # now look for consecutive intersections through the same edge, and remove them recursively
                # until all are gone (these correspond to bigons)
                def has_bigons(alg_ints_tuples):
                    for i in range(len(alg_ints_tuples)-1):
                        for j in range(i+1,len(alg_ints_tuples)):
                            if alg_ints_tuples[i][1] == alg_ints_tuples[j][1]:
                                return True, alg_ints_tuples, (i,j)
                    return False, alg_ints_tuples, (None,None)

                # A normal arc with endpoint at a vertex v of triangle t should have its first intersection
                # on the edge across from v. If the first intersection is with an edge that shares the vertex
                # v, then this gives a bigon (with one vertex v, hence and end bigon). If we see this we should
                # omit that intersection, such omission corresponds to isotoping across the bigon to remove it.
                def has_end_bigons(alg_ints_tuples):
                    if len(alg_ints_tuples) > 0:
                        if lifted_edge.vertex(0).index() in vertices[alg_ints_tuples[0][1]]:
                            return True, alg_ints_tuples, 0 
                        elif lifted_edge.vertex(1).point() in vertices[alg_ints_tuples[-1][1]]:
                            return True, alg_ints_tuples, -1
                    return False, alg_ints_tuples, None

                bigons, alg_ints_tuples, (i,j) = has_bigons(alg_ints_tuples)
                while bigons:
                    _ = alg_ints_tuples.pop(j)
                    _ = alg_ints_tuples.pop(i)
                    bigons, alg_ints_tuples, (i,j) = has_bigons(alg_ints_tuples)

                end_bigons, alg_ints_tuples, i = has_end_bigons(alg_ints_tuples)
                while end_bigons:
                    _ = alg_ints_tuples.pop(i)
                    end_bigons, alg_ints_tuples, i = has_bigons(alg_ints_tuples)

                for tup in alg_ints_tuples:
                    i = tup[1]
                    alg_arcs[i][lifted_edge.index()] += tup[0]
                    geom_arcs[i][lifted_edge.index()] += abs(tup[0])

            arcs = []
            for i in range(pc_tri.num_edges()):
                alg_arc = alg_arcs[i]
                geom_arc = geom_arcs[i]

                arc = Arc(self._combinatorial, vertices[i], alg_arc, geom_arc)
                arcs.append(arc)
            multi_arc = MultiArc(arcs)

        self._combinatorial._multi_arcs['pc_tri'] = multi_arc

        return self._combinatorial

#    def combinatorial_old(self):
#        if self._combinatorial_old is None:
#            triangles = []
#            edges_used = []
#            for e in self.edges():
#                for oriented_ind in [int(e.index()), ~int(e.index())]:
#                    if oriented_ind not in edges_used:
#                        e0 = e
#                        ind = e.index()
#                        triangle = [int(oriented_ind)]
#                        if oriented_ind >= 0:
#                            opp_vert = 1
#                        else:
#                            opp_vert = 0
#                        for i in range(2):
#                            v0 = e0.vertex(opp_vert)
#                            L = v0.valence()
#                            prev = (v0.incident_edges().index(ind)-1)%L
#                            ind = v0.incident_edges()[prev]
#                            e0 = self.edge(ind)
#                            if e0.vertex(0) == v0:
#                                triangle.append(int(ind))
#                                opp_vert = 1
#                            elif e0.vertex(1) == v0:
#                                triangle.append(~int(ind))
#                                opp_vert = 0
#                        v0 = e0.vertex(opp_vert)
#                        L = v0.valence()
#                        prev = (v0.incident_edges().index(ind)-1)%L
#                        ind = v0.incident_edges()[prev]
#                        assert ind == e.index()
#                        triangles.append(triangle)
#                        edges_used.extend(triangle)
#            vertex_markings = {}
#            for v in self.vertices():
#                if v.point() in self.top_poly().postcritical_set():
#                    vertex_markings[tuple(v.incident_edges())] = True
#                else:
#                    vertex_markings[tuple(v.incident_edges())] = False
#
#            self._combinatorial_old = CurverTriangulation(triangles,vertex_markings)
#        return self._combinatorial_old


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


    def plot(self, thickness=.7, DPI=400, aspect_ratio=1,show_vertices=True, arrow_size=.03, point_size=1, show_pc_tri=False, show_axes=False):
        G = Graphics()
        lines = []
        labels = []
        for i in range(self.num_edges()):
            edge = self.edge(i)
            points = edge.finite_points()

            edge_length = sum([seg.length() for seg in edge.segments()])
            seg_ind = 0
            L = edge.segment(seg_ind).length()
            while L < edge_length/2:
                seg_ind += 1
                L += edge.segment(seg_ind).length()
            seg = edge.segment(seg_ind)
            a = seg.endpoint(0).cif()
            b = seg.endpoint(1).cif()
            c = (a + b)/2
            a1 = (2*(b - c)/seg.length())*arrow_size*(exp(5*pi*I/6)) + c
            a2 = (2*(b - c)/seg.length())*arrow_size*(exp(-5*pi*I/6)) + c
            l = (2*(b - c)/seg.length())*1.5*arrow_size*(exp(pi*I/3)) + c



            ### for some reason Sage will not plot a line if its endpoints are elements of the ComplexIntervalField
            ### with imaginary part =0. So, as a somewhat hacky workaround we'll just add a very small multiple of I
            ### to each point.
            eps = CIF(0.00000000001*I)
            if self.is_lifted_tri():
                lines.append(Line([points[j].cif()+eps for j in range(len(points))], color=COLORS[edge.maps_to()], thickness=thickness, zorder=20))
                
                lines.append(Line([c,a1], color=COLORS[edge.maps_to()], thickness=thickness, zorder=20))
                lines.append(Line([c,a2], color=COLORS[edge.maps_to()], thickness=thickness, zorder=20))
                labels.append(Text(str(edge.index()),l, color=COLORS[edge.maps_to()], fontsize=6, zorder=20))                
                
                if show_pc_tri:
                    T = self.top_poly().pc_triangulation()
                    for i in range(T.num_edges()):
                        edge = T.edge(i)
                        points = edge.finite_points()
                        color = COLORS[i%12]
                        lines.append(Line([points[j].cif()+eps for j in range(len(points))], color=color, thickness=thickness*2, alpha=.5, zorder=15))
            else:
                f = self.num_finite_edges()
                color = COLORS[i%12]

                lines.append(Line([c,a1], color=color, thickness=thickness, zorder=20))
                lines.append(Line([c,a2], color=color, thickness=thickness, zorder=20))
                labels.append(Text(str(edge.index()),l, color=color, fontsize=6, zorder=20))                

                lines.append(Line([points[j].cif()+eps for j in range(len(points))], color=color, thickness=thickness, zorder=20))

        for L in lines:
            G += L
        for T in labels:
            G += T
        if show_vertices:
            for v in self.finite_vertices():
                color = 'black'
                for c in self.top_poly().postcritical_set():
                    if v.point().cif().overlaps(c.cif()):
                        color = 'red'
                        break
                p = (v.point().cif().real(), v.point().cif().imag())
                G += Disk(p, point_size*.03, (0,2*pi),  color='white', zorder=30)
                G += Disk(p, point_size*.03, (0,2*pi),  color=color, fill=False, thickness=thickness, zorder=35)
                G += Text(str(v.index()), p,  color=color, fontsize=5, zorder=35, vertical_alignment='center', horizontal_alignment='center')

        G.show(dpi=DPI, aspect_ratio=aspect_ratio, axes=show_axes)

class LiftedTriangulation(Triangulation):
    def __init__(self, marked_edges, top_poly):
        Triangulation.__init__(self, marked_edges, top_poly)

        self._pc_triangulation = top_poly.pc_triangulation()
        self._is_lifted = True
        self._vertices[-1] = MarkedVertex(Point(Infinity), self._pc_triangulation.vertices()[-1], self._vertices[-1].index())

#class CurverTriangulation(curverTriangulation):
#    def __init__(self,triangles,vertex_markings):
#        c_triangles = []
#        for t in triangles:
#            Ct = curverTriangle([curverEdge(t[i]) for i in range(len(t))])
#            c_triangles.append(Ct)
#
#
#        curverTriangulation.__init__(self,c_triangles)
#
#        self._vertex_markings = dict([(v,False) for v in self.vertices])
#        for v in self.vertices:
#            for w in vertex_markings:
#                if set(w) == set([e.index for e in v]):
#                    self._vertex_markings[v] = vertex_markings[w]
#
#    def vertex_markings(self):
#        return self._vertex_markings
#
#    def is_post_critical(self, vertex):
#        return self.vertex_markings()[vertex]

class Edge:
    def __init__(self, vertices, points, index = None):
        self._vertices = vertices
        self._points = points
        assert self._points[0] == self._vertices[0].point()
        assert self._points[-1] == self._vertices[1].point()

        if points[-1].alg() == Infinity:
            self._finite_points = points[:-1]
            self._ray = Ray([points[-2], points[-1]], self)
        else:
            self._finite_points = points
            self._ray = None
        self._segments = tuple([Segment([self._finite_points[i],self._finite_points[i+1]], self) for i in range(len(self._finite_points)-1)])
        self._index = index
        self._triangulation = None

    def triangulation(self):
        return self._triangulation

    def vert(self,i):
        return self._vertices[i]

    def vertex(self,i):
        return self._vertices[i]

    def vertices(self):
        return self._vertices

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

#    def alg_intersection(self,other):
#        # We define the algebraic intersection number as follows: for each transverse intersection 
#        # of self with other of the form shown, add +/- 1 as shown:
#
#        #                                        self
#        #                          ^    +1         |     -1
#        #                          |               |
#        #              other o-----|-------->------|---------o
#        #                          |               |
#        #                          |               v                   
#        #                        self
#        #
#        # If a vertex of self lies on the interior of other, then we add +/- 1/2 as shown:
#
#        #             self         self        other o------o-->---------o-------o
#        #               |           |                       |            |    
#        #               |           |                       |  +1/2      |  -1/2      
#        #               ^           v                       ^            v        
#        #          +1/2 |      -1/2 |                       |            |           
#        #               |           |                       |            |          
#        #   other o-----o-->--------o-----o                self         self                           
#        #
#        # In the case where self is an edge of the lifted triangulation, and other is an edge of the
#        # pc triangulation (which is the case we care about), it cannot happen that a vertex of other lies
#        # in the interior of self, so we don't need to worry about that case.
#
#        # we will assume that other is an edge that consists of a single segment, since that's all we
#        # need for now (all edges of the pc triangulation consist of a single segment).
#
#        assert other.num_segments() == 1
#
#        other_seg = other.segment(0)
#        alg_intersection = int(0)
#
#        segment_list = list(self.segments())
#
#        # if segment_list[0] has its initial point at a vertex of other, then we iteratively remove segments 
#        # from the beginning of segment_list until segment_list[0] has its initial point disjoint from other.
#        # Similarly if segment_list[-1] has its terminal point at a vertex of other, then we iteratively remove
#        # segments from the end of segment_list until segment_list[-1] has its terminal vertex disjoint from
#        # other.
#
#        # Once we have adjusted segment_list as described above, we can determine the algebraic intersection
#        # number of self with other by adding up intersection numbers of individual segments of self. 
#
#        if segment_list[0].endpoint(0) == other.point(0) or segment_list[0].endpoint(0) == other.point(1):
#            while other_seg.on_seg(segment_list[0].endpoint(0)):
#                _ = segment_list.pop(0)
#                if len(segment_list) == 0:
#                    return 0
#
#        if segment_list[-1].endpoint(1) == other.point(0) or segment_list[-1].endpoint(1) == other.point(1):
#            while other_seg.on_seg(segment_list[-1].endpoint(1)):
#                _ = segment_list.pop(-1)
#                if len(segment_list) == 0:
#                    return 0
#
#
#
#        for seg in segment_list:
#            if other_seg.on_seg(seg.endpoint(0)):
#                if other_seg.below_seg(seg.endpoint(1)):
#                    alg_intersection -= 1
#                elif other_seg.above_seg(seg.endpoint(1)):
#                    alg_intersection += 1
#            elif other_seg.below_seg(seg.endpoint(0)):
#                if other_seg.on_seg(seg.endpoint(1)):
#                    alg_intersection += 1
#                elif other_seg.above_seg(seg.endpoint(1)):
#                    if other_seg.transverse_to(seg):
#                        alg_intersection += 2
#            elif other_seg.above_seg(seg.endpoint(0)):
#                if other_seg.on_seg(seg.endpoint(1)):
#                    alg_intersection -= 1
#                elif other_seg.below_seg(seg.endpoint(1)):
#                    if other_seg.transverse_to(seg):
#                        alg_intersection -= 2
#
#        geom_int = abs(alg_intersection/2)
#        assert geom_int == int(geom_int)
#
#        return int(geom_int)
#
#
#
#

class MarkedEdge(Edge):
    def __init__(self, vertices, points, maps_to, index=None):
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

    def slope(self): ## return an RIF, not an algebraic integer
        a,b = self.endpoint(0).cif(), self.endpoint(1).cif()
        if self._slope == None:
            if a.real().overlaps(b.real()):
                self._slope = Infinity
            else:
                self._slope = (a.imag() - b.imag())/(a.real() - b.real())
        return self._slope

    def unit_tangent(self):
        tangent = self.endpoint(1).alg() - self.endpoint(0).alg()
        unit_tangent = tangent/tangent.abs()
        return unit_tangent

    def y_intercept(self): ## return an RIF, not an algebraic integer
        if self.slope() != Infinity:
            return -self.slope()*self.endpoint(0).cif().real() + self.endpoint(0).cif().imag()
        else:
            return None

    def transverse_to(self, other):
        '''Return True if self has nontrivial transverse intersection with other.
        '''
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
        '''Return True if point lies on this segment (note it must lie on the segment, not just
            the line that segment is contained in.)
        '''
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

    def length(self):
        return abs(self.endpoint(0).cif() - self.endpoint(1).cif())


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

    def degree(self):
        return len(self.incident_edges())

    def incident_edges(self):
        return self._incident_edges

    def __eq__(self,other):
        return self.point() == other.point()

    def __ne__(self,other):
        return self.point() != other.point()

    def __repr__(self):
        return self.point().__repr__()

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

    def _set_alg(self, new_alg):
        self._alg = new_alg
        self._set_cif()


    def cif(self):
        return self._cif

    def _set_cif(self):
        if self._alg == Infinity:
            self._cif = CIF(Infinity)
        else:
            self._cif = self._alg.interval(ComplexIntervalField(self._precision))

    def precision(self):
        return self._precision

    def set_precision(self, prec):
        self._precision = prec
        self._cif = self._alg.interval(ComplexIntervalField(prec))

    def overlaps(self,other):
        return self.cif().overlaps(other.cif())

    def __eq__(self,other):
        return type(self) == type(other) and self.overlaps(other)

    def __ne__(self,other):
        return not self.__eq__(other) 

    def __repr__(self):
        return self.cif().__repr__()

    def __str__(self):
        return self.cif().__str__()

    def __mul__(self,other):
        return Point(self.alg()*other)

    def __add__(self,other):
        return Point(self.alg() + other)

    def __sub__(self,other):
        return Point(self.alg() - other)





