from __future__ import print_function
from lifty.sage_types import *
import random
import itertools as it
import sys
import copy

from lifty.kernel.triangulation import Point, Segment, Edge, Vertex, Triangulation, LiftedTriangulation
from lifty.constants import RABBIT, CORABBIT, AIRPLANE
from lifty.kernel.errors import IntervalError, IterationError, SubdivisionError
from lifty.kernel.combinatorial_tri import Arc, MultiArc

COLORS = ['blue','red','green','cyan','black','grey','pink','purple','orange','brown','limegreen','olive']

class TopPoly:
    def __init__(self, coefficients):
        """ A Polynomial is specified by a list (or a tuple) of coefficients [a_n, a_{n-1}, ..., a_0]. 
            Each of which is an algebraic number. For example, to enter the rabbit:
            
            sage: x = polygen(ZZ)
            sage: c = QQbar.polynomial_root((x^2+x)^2+x, CIF(RIF(-0.13,-0.12),RIF(0.74,0.75)))
            sage: P = TopPoly([1,0,c])

        """
        self._coefficients = coefficients

        self._max_precision = 1000

        self.degree = len(self._coefficients)-1
        self._roots = None
        self._critical_pts = None
        self._postcritical_set = None
        self._pc_triangulation = None
        self._lifted_pcs = None
        self._lifted_tri = None
        self._derivative = None

    def max_precision(self):
        return self._max_precision

    def set_max_precision(self, bits):
        self._max_precision = bits

    def coefficients(self):
        return self._coefficients
    
    def coeff(self, k):
        return self._coefficients[self.degree-k]

    def evaluate(self,z):
        result = 0
        for i in range(self.degree+1):
            result += self.coeff(i)*(z**i)
        return result

    def iterate(self,z,k):
        """ evaluate P^k(z)
        """
        result = self.evaluate(z)
        for i in range(k-1):
            result = self.evaluate(result)
        return result

    def derivative(self):
        """ return the derivative of self, as a TopPoly.
        """
        if self._derivative == None:
            Q_coeffs = []
            P_coeffs = self.coefficients()
            for i in range(self.degree):
                n = self.degree - i
                Q_coeffs.append(n*P_coeffs[i])
            self._derivative = TopPoly(Q_coeffs)
        return self._derivative

    def sage_poly(self):
        """ return a Sage polynomial with coefficients self.coefficients
        """
        d = self.degree
        x = polygen(QQbar)
        p=self.coeff(d)*x**d
        for m in range(d):
            p += self.coeff(m)*x**m
        return p

    def _find_roots(self):
        return complex_roots(self.sage_poly(),retval='algebraic')

    def roots(self):
        """ return the roots of self. Results are cached.
        """
        if self._roots == None:
            self._roots = self._find_roots()
        return self._roots

    def critical_pts(self):
        """ return the critical points of self. Results are cached.
        """
        if self._critical_pts == None:
            Q = self.derivative()
            Qroots = Q.roots()
            self._critical_pts = [Point(Qroots[i][0]) for i in range(len(Qroots))]
        return self._critical_pts

    def postcritical_set(self, max_iterations=20):
        """ return the postcritical set of self. Results are cached. 

            sage: P=lifty.quadratic('rabbit')
            sage: P.postcritical_set(max_iterations=20)

                [-0.6623589786223730? + 0.5622795120623013?*I,
                 -0.1225611668766537? + 0.744861766619745?*I,
                 0]

            The argument max_iterations defaults to 20. If a critical point has period > 20, this will need to be increased.
        """
        if self._postcritical_set == None:
            M = max_iterations
            C = [c.alg() for c in self.critical_pts()]
            PCS = dict()
            for c in C:
                confirmed_finite = False
                PCS[c] = []
                pcp = self.evaluate(c)
                PCS[c].append(pcp)
                for i in range(1,M+1):
                    pcp = self.evaluate(pcp)
                    if pcp not in PCS[c]:
                        PCS[c].append(pcp)
                        continue
                    else:
                        confirmed_finite = True
                        PCS[c].reverse()  ## ensures that the critical point, if in PCS[c], is the first in the list
                        break
                if confirmed_finite == False:
                    raise IterationError('either this polynomial is not post-critally finite, or max_iterations needs to be increased.')
                else:
                    continue
            postcritical = list(set([pcp for c in C for pcp in PCS[c]]))
            pcs = [Point(p) for p in postcritical]

            pcs.sort(key = lambda x: x.cif().real())

            # check that the intervals for points in the postcritical set are all disjoint. If not, increase precision.
            make_intervals_disjoint(pcs, self.max_precision())
            self._postcritical_set = tuple(pcs)
        return self._postcritical_set



    def pc_triangulation(self):
        """Return a triangulation whose vertices are the postcritical set. This will always be a piecewise
            linear triangulation (each edge is a single line segment). Note that at this point 
            we don't include edges with endpoints at infinity. Results are cached.
        """
        if self._pc_triangulation == None:
            pcs = self.postcritical_set()
            vertices = [Vertex(p) for p in pcs]
            L = len(pcs)
            edge_list = []
            for i in range(L-1):
                for j in range(i+1,L):
                    a,b = pcs[i], pcs[j]
                    v,w = vertices[i], vertices[j]
                    e = Edge([v,w],[a,b])
                    s = e.segment(0)
                    if e not in edge_list:
                        to_trash = False
                        for edge in edge_list:
                            if s.transverse_to(edge.segment(0)):
                                to_trash = True
                                break
                        if not to_trash:
                            for k in range(L):
                                if k != i and k != j:
                                    if s.on_seg(pcs[k]):
                                        to_trash = True
                                        break
                        if not to_trash:
                            edge_list.append(e)

            for i in range(len(edge_list)):
                edge_list[i]._index = i

            T = Triangulation(edge_list, self)
            T.add_edges_to_infinity()
            self._pc_triangulation = T

        return self._pc_triangulation




    def lifts(self,point):
        """Return all lifts of point, which must be a Point.
        """
        c = self.coeff(0)-point.alg()
        coeffs = self.coefficients()[:-1]+[c]
        Q = TopPoly(coeffs)
        Qroots = Q.roots()
        lifts = [Point(Qroots[i][0]) for i in range(len(Qroots))]
        return lifts

    def lifted_pcs(self):
        """Return lifts of the postcritical set. Results are cached.
        """
        if self._lifted_pcs == None:
            lifted_pcs = []
            for p in self.postcritical_set():
                lifts = self.lifts(p)
                for l in lifts:
                    lifted_pcs.append(l)
            lpcs = lifted_pcs
            make_intervals_disjoint(lpcs, self.max_precision())

            self._lifted_pcs = lpcs
        return self._lifted_pcs


    def dist_to_pcs(self,segment):
        """Returns minimal distance between the given edge and points of the 
            postcritical set which are not vertices of the edge. This should probably
            be moved out of the TopPoly class eventually.
        """
        a,b = segment[0].cif(), segment[1].cif()
        pcs = [p for p in self.postcritical_set() if p not in segment]
        m = Infinity
        for i in range(len(pcs)):
            p = pcs[i].cif()
            x_0,y_0 = p.real(), p.imag()
            b_re, b_im = b.real(), b.imag()
            a_re, a_im = a.real(), a.imag()
            dx = b_re - a_re
            dy = b_im - a_im
            c = -2*(x_0*dx+y_0*dy) + 2*(dx*a_re+dy*a_im)
            s = 2*(dx**2 + dy**2)
            t = -c/s
            if t.overlaps(RIF([0,1])):
                d = abs(p - b*t + (1-t)*a)
            else:
                d = min([abs(p - a), abs(p - b)])
            if d < m:
                m = d
            elif d.overlaps(m):
                m = m.union(d)
            else:
                continue
        return m

## Let p(z) be a complex polynomial. Then by the mean value thm, Re( (p(B)-p(A))/(B-A) ) <= max(Re(p'(z))),
## and Im( (p(B)-p(A))/(B-A) ) <= max(Im(p'(z))), where the maximum is taken over the segment 
## from A to B. This implies that | ((p(B) - p(A))/(B - A) | <= 2*max(|p'(z)|).

## Now let a,b be points in the complex plane, and let A and B be lifts of a and b by p. 
## Denote by a-->b the line segment from a to be, and by A-->B the segment from A to B. We
## want to check that some lift of a-->b is isotopic to A-->B. Or, equivalently, that
## p(A-->B) is isotopic to a-->b. If this fails, then p(A-->B) must be an arc that has 
## length at least twice d, where d is the distance from a-->b to {postcritical set}\{a,b}.
## It follows that there is some point C on A-->B such that the length of a-->p(C) => d,
## i.e., |P(A)-P(C)| >= d. Thus we get |B-A| >= |C-A| >= d/(max(|p'(z)|)), where the maximum is taken
## over the segment from A to B. The quantity d is computed by the above function, dist_to_pcs().

## Let p'(z) = a_n*z^n+a_{n-1}*z^{n-1}+ ... + a_1*z + a_0. Then 
## max(|p'(z)|) <= max( |a_n*z^n| + |a_{n-1}*z^{n-1}| + ... + |a_1*z| + |a_0| ). The right hand side
## is maximized at the point on A-->B with largest modulus, which must be one of the endpoints since it 
## is a line segment. 

## Note that we do all of the above to avoid having to maximize a polynomial, at the cost of most likely 
## getting a worse bound on | ((p(B) - p(A))/(B - A) | than could have been obtained otherwise. My guess 
## is that this is worth it, but maybe there are improvements to be made here.



    def max_deriv(self,segment):
        """If Q(x) = a_n*z^n+a_{n-1}*z^{n-1}+ ... + a_1*z + a_0 is the derivative of self,
            then return the maximum of |a_n*z^n|+|a_{n-1}*z^{n-1}|+ ... + |a_1*z| + |a_0|
            over the segment, which will be attained at an endpoint.  
            This should probably be moved out of the TopPoly class eventually.
        """
        Q = self.derivative()
        s = max([abs(segment[0].cif()), abs(segment[1].cif())])

        d = Q.degree
        max_deriv = 0
        for i in range(d+1):
            a = Q.coeff(i)
            max_deriv += abs(a)*s**i
        return max_deriv


    def lift_edge(self,edge):
        """Return all lifts of edge, as a list of lists of algebraic numbers.
        """
        Q = self.derivative()
        a,b = edge.point(0), edge.point(1)
        a_lifts, b_lifts = self.lifts(a), self.lifts(b)
        make_intervals_disjoint(a_lifts, self.max_precision())
        make_intervals_disjoint(b_lifts, self.max_precision())
        
        lifted_edges = []

        ## to determine how many edges should have their initial point at each of the lifts of a,
        ## we lift a point c along the line from a to b that is very close to a, so that each lift of 
        ## c is distance <D from some lift of a. For a given lift A of a, the number of edges emanating 
        ## from A is the number of lifts of c that are D-close to A.
        t = QQbar(1)/100000
        def init_edge(t, lifted_edges):
            c = Point(b.alg()*t + (1-t)*a.alg())
            c_lifts = self.lifts(c)
            for A in a_lifts:
                for C in c_lifts:
                    D = self.dist_to_pcs([a,c])/(2*self.max_deriv([A,C]))
                    if abs(A.cif()-C.cif()) < D:
                        lifted_edges.append([A])
            return t, lifted_edges

        # if t is too large, we may not get all edges, so we keep making t smaller until we get self.degree edges
        while len(lifted_edges) < self.degree:
            t, lifted_edges = init_edge(t/10, lifted_edges)

        assert len(lifted_edges) == self.degree # this should always be true, but could fail if there is a bug or bad math.

        ## we extend the lift of each edge one segment at a time, until we have the full lift. Each
        ## iteration of below function extends by one segment. 
        def extend_edges(t, dt, interval, lifted_edges, finished):
            s1 = Point(b.alg()*t + (1-t)*a.alg())
            s2 = Point(b.alg()*(t+dt)+(1-(t+dt))*a.alg())
            lifted_segs = []
    
            seg = [s1,s2]
            s1_lifts = []
            for e in lifted_edges:
                if e[-1] not in s1_lifts:
                    s1_lifts.append(e[-1])
            s2_lifts = self.lifts(s2)
            make_intervals_disjoint(s2_lifts, self.max_precision())
    
            for S1,S2 in it.product(s1_lifts,s2_lifts):
                
                ## if A and B are lifts of a and b, and abs(A-B) < D, then the image under self of the segment A--->B 
                ## is isotopic to the segment a--->b. So the segment A--->B is isotopic to a lift of a--->b.
                D = self.dist_to_pcs(seg)/(2*self.max_deriv([S1,S2]))
                if abs(S1.cif()-S2.cif()) < D:
                    lifted_segs.append([S1,S2])
            
            ## if we have fewer than self.degree lifted segments, then we are missing some. This would
            ## be caused by trying to lift a subarc of the edge that is too long. So we decrease dt, and 
            ## adjust the interval inside which we know the correct dt should live.
            if len(lifted_segs) < self.degree:
                interval = [interval[0], dt]
                dt = (interval[0]+dt)/2
                return t, dt, interval, lifted_edges, finished

            ## if we have more than self.degree lifted segments, then we must be close to a critical point, and 
            ## we are lifting a subarc of the edge that is too short (Since we are near a critical point, multiple
            ## lifted subarcs are close together. By increasing the length of the lift, we can ensure that an arc
            ## between a lift of the initial point and a lift of the terminal point is actually a lift of the subarc). 
            ## So we increase dt, and adjust the interval inside which we know the correct dt should live.
            elif len(lifted_segs) > self.degree:
                interval = [dt, interval[1]]
                dt = (dt+interval[1])/2
                return t, dt, interval, lifted_edges, finished

            ## if we get self.degree lifts, then we are happy. Set the new t for the next iteration to be t+dt.
            ## We set dt to something that (we hope) is a good guess for what should work. It shouldn't change much
            ## from what worked for the previous segment, but making it 2*dt (instead of dt) seems to work pretty well.
            else:
                t = t + dt
                dt = 2*dt
                interval = [0,1-t]

                for i in range(len(lifted_edges)):
                    lifted = lifted_edges[i]
                    if lifted[-1] not in b_lifts:
                        if len(lifted) > 1: # without this (or some other solution) an edge can lift to multiple copies of one of its lift (if edge.vert(0) lift to a critical point)
                            for B in b_lifts:
                                D = self.dist_to_pcs([s1,b])/(2*self.max_deriv([lifted[-1],B]))
                                
                                ## try to jump to the end. If the segment from the last point in the lifted
                                ## edge so far to one of the lifts of b has length less than D, then this segment
                                ## will do, and we can finish.
                                if abs(lifted[-1].cif() - B.cif()) < D:
                                    lifted.append(B)
                                    finished[i] = True
                                    break
                        ## if we are unable to jump to the end, find the lifted segment that should be joined to 
                        ## the edge in progress to extend it.            
                        if not finished[i]:
                            for i in range(len(lifted_segs)):
                                seg = lifted_segs[i]
                                if lifted[-1].cif().overlaps(seg[0].cif()):
                                    lifted.append(seg[1])
                                    _ = lifted_segs.pop(i)
                                    break
            return t, dt, interval, lifted_edges, finished

        ## the first iteration of function for extend edges
        t, dt, interval, lifted_edges, finished = extend_edges(0, QQ(1)/10, [0,1], lifted_edges, [False for _ in range(self.degree)])
        print('|',end='')

        ## continue to iterate the extend_edge function until all edges are finished.
        while False in finished:
            t, dt, interval, lifted_edges, finished = extend_edges(t, dt, interval, lifted_edges, finished)
            print('|',end='')
        return lifted_edges

    def lifted_tri(self):
        return self.lifted_triangulation()

    def lifted_triangulation(self):
        """Return the lift of the triangulation self.pc_triangulation(). As self.pc_triangulation() does not include 
            edges with endpoint at infinity, this one won't either. Results are cached.
        """
        if self._lifted_tri == None:
            T = self.pc_triangulation()
            all_edges = []

            # lift all the edges
            vertices = []
            subvertices = []

            # make vertices. These are MarkedVertices, so that we store info about where each one maps to in the pc_triangulation
            pc_verts = [v for v in T.finite_vertices()]
            for v in pc_verts:
                p = v.point()
                p_lifts = self.lifts(p)
                for P in p_lifts:
                    vertices.append(Vertex(P,v.index()))


            vertex_pts = [v.point() for v in vertices]
            for i in range(T.num_finite_edges()):
                e = T.edge(i)
                e_lifted_edges = self.lift_edge(e)
                print('edge {} lifted.'.format(e.index()))

                for points in e_lifted_edges:
                    j,k = vertex_pts.index(points[0]), vertex_pts.index(points[-1])
                    v,w = vertices[j], vertices[k]
                    points[0] = vertices[j].point()
                    points[-1] = vertices[k].point()
                    subvertices += points[1:-1]
                    edge = Edge([v,w],points, i, None)
                    all_edges.append(edge)


            # set indices so they agree with the order in which edges are listed.
            for i in range(len(all_edges)):
                e = all_edges[i]
                e._index = i

            lpcs = self.lifted_pcs()

            # precision required for points in the lifted postcritical set to have disjoint intervals.
            prec = max([p.precision() for p in lpcs]) 

            # increase the precision for vertices to prec, so we can compare intervals to get rid of duplicates.
            for v in vertex_pts:
                v.set_precision(prec)


            # now we can make sure all points in the triangulation have disjoint intervals.
            make_intervals_disjoint(subvertices+vertex_pts, self.max_precision())

            # make sure the length of the shortest segment is large compared to the precision
            short_seg = min([s.length() for e in all_edges for s in e.segments()])
            min_prec = min([p.precision() for p in subvertices+vertex_pts])
            dec_prec = 10**(-(min_prec*log(2,10)))
            if short_seg < dec_prec*1000000:
                new_prec = int((-log(short_seg,10)+7)/log(2,10))
                for p in subvertices+vertex_pts:
                    if p.precision() < new_prec:
                        p.set_precision(new_prec)


            T = LiftedTriangulation(all_edges, self)

            is_embedded, bad_segs = T.is_embedded()
            while not is_embedded:
                for bad_seg in bad_segs:
                    edge_index = bad_seg[0]
                    seg_index = bad_seg[1]
                    edge = T.edge(edge_index)
                    seg = edge.segment(seg_index)
                    S0 = seg[0]
                    S1 = seg[1]
                    s0, s1 = Point(self.evaluate(S0.alg())), Point(self.evaluate(S1.alg()))
                    c = (s0.alg() + s1.alg())/2
                    c_lifts = self.lifts(c)
                    D = self.dist_to_pcs([s0,s1])/(2*self.max_deriv(seg))
                    candidates = [C for C in c_lifts if abs(C.cif() - S0.cif()) < D]
                    if len(candidates) == 1:
                        C = candidates[0]
                        T.subdivide_segment(edge_index, seg_index, C)
                    else:  ## I don't think this should ever happen, but just to be safe...
                        raise SubdivisionError('unable subdivide edge.')
                is_embedded, bad_segs = T.is_embedded()
            T.add_edges_to_infinity()
            self._lifted_tri = T
        return self._lifted_tri

    def invariant_tree(self):
        

        T = self.pc_triangulation()
        L = self.lifted_triangulation()
        TCom = T.combinatorial()
        LCom = L.combinatorial()

        LCom_collapsed, collapse_list = LCom.collapse_to_ideal()

        LCom_flipped_to_pc, flip_list = LCom_collapsed.flip_to_pc_tri()

        # LCom_flipped_to_pc is isometric to T, but with edges labelled differently.
        # So we need a dictionary to maps each edge of LCom_flipped_to_pc to the isotopic 
        # edge of T. We also need the dictionary to tell us if the two edges are oriented the same.
        # We use the fact that lifting, collapsing and flipping are done in such a way that the 
        # vertex indices of post-critical vertices is preserved. So the indexing of vertices 
        # will be the same for T, L, TCom, LCom, LCom_collapsed, and LCom_flipped_to_pc.
        # Since T is constructed is such a way that edges of T are uniquely defined by their
        # vertices, we can use this to identify edges of T with edges of LCom_flipped_to_pc.
        
        edge_map = {}

        T_edge_verts = {(e.vertex(0).index(), e.vertex(1).index()):e for e in T.edges()}
        LCom_flipped_edge_verts = {(e.vertex(0).index(), e.vertex(1).index()):e for e in LCom_flipped_to_pc.edges()}

        for (a,b) in LCom_flipped_edge_verts:
            if (a,b) in T_edge_verts:
                edge_map[LCom_flipped_edge_verts[(a,b)].index()] = (T_edge_verts[(a,b)].index(),True)
            else:
                edge_map[LCom_flipped_edge_verts[(a,b)].index()] = (T_edge_verts[(b,a)].index(),False)

        def transfer_multi_arc(CTri, edge_map, multi_arc):
            intersection_sequence = []
            M = multi_arc
            arcs = []
            for arc in M.arcs():
                intersection_sequence = []
                for e in arc.intersection_sequence():
                    sign = e.sign()
                    if edge_map[e.index()][1] == False:
                        sign *= -1
                    if sign == 1:
                        intersection_sequence.append(edge_map[e.index()][0])
                    elif sign == -1:
                        intersection_sequence.append(-edge_map[e.index()][0]-1)
                arcs.append(Arc(CTri,[arc.vertex(0).index(),arc.vertex(1).index()],[],[],intersection_sequence))
            M0 = MultiArc(arcs)
            CTri.add_multi_arc(M0,'tree')


        def iterate_lifting(TCom, LCom, collapse_list, flip_list, edge_map):

            # we will work with a copy of LCom
            L = copy.deepcopy(LCom)

            #first we lift the multi_arc 'tree' to LCom
            M = TCom.multi_arc('tree')
            arcs = []
            for arc in M.arcs():
                algebraic = arc.algebraic()
                geometric = arc.geometric()
                lifted_algebraic = [0 for i in range(L.num_edges())]
                lifted_geometric = [0 for i in range(L.num_edges())]
                for i in range(L.num_edges()):
                    e = L.edge(i).maps_to()
                    lifted_geometric[i] = geometric[e]
                    lifted_algebraic[i] = algebraic[e]
                arcs.append(Arc(L,[L.num_vertices()-1, L.num_vertices()-1], lifted_algebraic, lifted_geometric))
            M_lifted = MultiArc(arcs)
            L.add_multi_arc(M_lifted,'tree')

            print(L)
            print(L.multi_arc('tree'))

            # now collapse LCom to an ideal triangulation
            for edge in collapse_list:
                L.collapse_edge(edge)

            # now flip LCom to pc_triangulation
            for edge in flip_list:
                L.flip_edge(edge)

            # now transfer the multi_arc 'tree' from L to TCom.
            transfer_multi_arc(TCom, edge_map, L.multi_arc('tree'))

            return TCom

        a0 = Arc(TCom,[TCom.num_vertices()-1,TCom.num_vertices()-1],[],[],[-1,-3,-2])
        a1 = Arc(TCom,[TCom.num_vertices()-1,TCom.num_vertices()-1],[],[],[1,-5])
        a2 = Arc(TCom,[TCom.num_vertices()-1,TCom.num_vertices()-1],[],[],[4,2,3])
        a3 = Arc(TCom,[TCom.num_vertices()-1,TCom.num_vertices()-1],[],[],[-4,0])
        a4 = Arc(TCom,[TCom.num_vertices()-1,TCom.num_vertices()-1],[],[],[-4,-3,-2])
        M = MultiArc([a0,a1,a2,a3,a4])
        TCom.add_multi_arc(M,'tree')
        print(TCom)
        print(TCom.multi_arc('tree'))

        for i in range(5):
            TCom = iterate_lifting(TCom, LCom, collapse_list, flip_list, edge_map)
            print(TCom)
            print(TCom.multi_arc('tree'))
        #def relabel_edges(CTri, edge_map):
        #    # fix _signs for each triangle according to orientation of edges of T
        #    for triangle in CTri.triangles():
        #        signs = list(triangle._signs)
        #        for i in range(3):
        #            if not edge_map[triangle.edge(i).index()][1]:
        #                signs[i] = triangle._signs[i]*(-1)
        #        triangle._signs = tuple(signs)
        #    # clear incident edges for each vertex (will be recomputed on next call)
        #    for v in CTri.vertices():
        #        v._incident_edges = None
    #
        #    # re-index edges of LCom_flipped_to_pc according to edge_map
        #    for edge in CTri.edges():
        #        if not edge_map[edge.index()][1]:
        #            v0,v1 = edge.vertex(0), edge.vertex(1)
        #            edge._vertices = (v1,v0)
        #        #clear adjacent triangles of positive and negative reps (will be recomputed on next call)
        #        edge.positive_rep()._adjacent_triangle = None
        #        edge.negative_rep()._adjacent_triangle = None
        #        edge._index = edge_map[edge.index()][0]
#
        #    # fix multi-arc edge intersection info
        #    num_edges = CTri.num_edges()
        #    for key in CTri.multi_arcs():
        #        M = CTri.multi_arc(key)
        #        for a in M.arcs():
        #            intersection_sequence = []
        #            for e in a.intersection_sequence():
        #                if edge_map[e.index()][1] == True:
        #                    intersection_sequence.append(e)
        #                else:
        #                    intersection_sequence.append(e.opp_signed_edge())
        #            a._intersection_sequence = intersection_sequence
        #            a._recompute_intersection_vectors()


        return TCom


        # Now LCom_flipped_to_pc should have edges with indices and orientations the same as T

        # We need to create a filling multi-arc with all vertices at infinity in TCom, to initialize
        # the invariant tree. There is an easy way to do this: For each vertex v of TCom with an edge to 
        # infinity (such a vertex is distance 1 from inf), put an arc from inf, around v, and back to inf. 
        #Then for each vertex v that has an edge to a distance 1 vertex (such a vertex is distance 2 from inf),
        # put an edge that tracks first the edge from infinity, then the edge to v, then goes around v 
        # and tracks the two edges back to infinity. Continue in this way (i.e., next do all vertices 
        # that are distance 3 from infinity, etc.). When this process terminates, some subset of the arcs 
        # are "outermost arcs" with respect to this iterative process, and these outermost arcs are cyclically
        # ordered around the vertex at infinity, and we may number them 1 to n based on this ordering (note
        # that this numbering is only well-defined up to cyclic permutation). Now for each even k, 1<k<n, 
        # put a new arc starting between arc n and arc 1, and ending between arc k and k+1. For a given k,
        # such an arc will intersect every edge that has one vertex among those enclosed by arcs 1,2,...,k,
        # and another vertex among those enclosed by arcs k+1,...,n.




            




        ###
        # 1. compute the pc_triangulation and lifted triangulation, then combinatorilize
        # 2. Find sequence of collapses and flips that take the lifted_tri to the pc_tri
        # 3. Encode the lifting map on arcs combinatorially
        # 4. choose a filling arc system A and encode it in combinatorial pc_tri.
        # 5. Lift A to lifted_tri, track it through collapses/flips to get A'
        # 6. Check to see if the Hubbard tree is in a 2-nbhd of A'.
        # 7. If #6 gives No, check for obstruction.
        # 8. If #7 gives No, rename A' to A, go to #5.


def make_intervals_disjoint(points, max_precision):
    for i in range(len(points)-1):
        for j in range(i+1,len(points)):
            p, q = points[i], points[j]
            prec = max([p.precision(), q.precision()])
            while p.overlaps(q):
                prec = 2*prec
                if prec < max_precision:
                    p.set_precision(prec)
                    q.set_precision(prec)
                else:
                    raise IntervalError('max precision reached: two of the points are too close together to be distinguished by intervals.')


def get_random_tp(degree,period):
    """Return a random polynomial P(z) = z^degree + c, for which the critical point 0 is (period)-periodic
        (i.e., P^{period}(0)=0). There are degree^{period-1} of these, so the returned value is random in the
        sense that it chooses one of these randomly.
    """
    x = polygen(QQbar)
    d = degree
    p = x**d + x
    for i in range(period-2):
        p = p**d + x
    roots = [r[0] for r in complex_roots(p,retval='algebraic')]
    roots.remove(0)
    c = random.choice(roots)
    coeffs = [1]+[0 for i in range(d-1)]+[c]
    return TopPoly(coeffs)

def get_example(degree,period,root):
    x = polygen(QQbar)
    d = degree
    p = x**d + x
    for i in range(period-2):
        p = p**d + x
    roots = [r[0] for r in complex_roots(p,retval='algebraic')]
    roots.remove(0)
    c = roots[root%(len(roots))]
    coeffs = [1]+[0 for i in range(d-1)]+[c]
    return TopPoly(coeffs)


def quadratic(type):
    """ Return one of the three 3-periodic quadratic polynomials z^2+c. The argument "type" can be
        one of 'rabbit', 'corabbit', or 'airplane'.
    """
    if type == 'rabbit':
        return TopPoly([1,0,RABBIT])
    elif type == 'corabbit':
        return TopPoly([1,0,CORABBIT])
    elif type == 'airplane':
        return TopPoly([1,0,AIRPLANE])
    else:
        return 'Error: type must be \'rabbit\', \'corabbit\', or \'airplane\'.'



#def relabel_edges(CTri, edge_map):
#    # fix _signs for each triangle according to orientation of edges of T
#    for triangle in CTri.triangles():
#        edges = triangle.edges()
#        signs = list(triangle._signs)
#        for i in range(3):
#            if not edge_map[triangle.edge(i).index()][1]:
#                signs[i] = triangle._signs[i]*(-1)
#        triangle._signs = tuple(signs)
#        triangle._initialize_signed_edges_and_vertices(edges)
#    # clear incident edges for each vertex (will be recomputed on next call)
#    for v in CTri.vertices():
#        v._incident_edges = None
#    # re-index edges of LCom_flipped_to_pc according to edge_map
#    for edge in CTri.edges():
#        if not edge_map[edge.index()][1]:
#            v0,v1 = edge.vertex(0), edge.vertex(1)
#            edge._vertices = (v1,v0)
#        #clear adjacent triangles of positive and negative reps (will be recomputed on next call)
#        edge.positive_rep()._adjacent_triangle = None
#        edge.negative_rep()._adjacent_triangle = None
#        edge._index = edge_map[edge.index()][0]
#    # fix multi-arc edge intersection info
#    num_edges = CTri.num_edges()
#    for key in CTri.multi_arcs():
#        M = CTri.multi_arc(key)
#        for a in M.arcs():
#            intersection_sequence = []
#            for e in a.intersection_sequence():
#                if edge_map[e.index()][1] == True:
#                    intersection_sequence.append(e)
#                else:
#                    intersection_sequence.append(e.opp_signed_edge())
#            a._intersection_sequence = intersection_sequence
#            a._recompute_intersection_vectors()




