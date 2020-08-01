
from lifty.sage_types import *
import random
import itertools as it
import sys

from lifty.kernel.triangulation import Point, Segment, Edge, MarkedEdge, Vertex, Triangulation, LiftedTriangulation
from lifty.constants import RABBIT, CORABBIT, AIRPLANE

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


            # check that the intervals for points in the postcritical set are all disjoint. If not, increase precision.
            make_intervals_disjoint(pcs, self.max_precision())
            self._postcritical_set = pcs
        return self._postcritical_set


    def pc_triangulation(self):
        """Return a triangulation whose vertices are the postcritical set. This will always be a 
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
                    e = Edge([v,w],[a,b], self)
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
        print(n(t),n(dt))

        ## continue to iterate the extend_edge function until all edges are finished.
        while False in finished:
            t, dt, interval, lifted_edges, finished = extend_edges(t, dt, interval, lifted_edges, finished)
            print(n(t),n(dt))
        return lifted_edges

    def lifted_tri(self):
        """Return the lift of the triangulation self.pc_triangulation(). As self.pc_triangulation() does not include 
            edges with endpoint at infinity, this one won't either. Results are cached.
        """
        if self._lifted_tri == None:
            T = self.pc_triangulation()
            all_edges = []

            # lift all the edges
            subvertices = []
            vertex_pts = self.lifted_pcs()
            vertices = [Vertex(p) for p in vertex_pts]
            for i in range(T.num_finite_edges()):
                e = T.edge(i)
                e_lifted_edges = self.lift_edge(e)
                print('edge {} lifted.'.format(e.index()))

                for points in e_lifted_edges:
                    j,k = vertex_pts.index(points[0]), vertex_pts.index(points[-1])
                    v,w = vertices[j], vertices[k]
                    subvertices += points[1:-1]
                    edge = MarkedEdge([v,w],points, i, None)
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

            self._lifted_tri = T
        return self._lifted_tri

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
                    print(p,q)
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
    roots = complex_roots(p,retval='algebraic')
    c = random.choice(roots)[0]
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

class SubdivisionError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'SubdivisionError: {0}'.format(self.message)

class IntervalError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'IntervalError: {0}'.format(self.message)

class IterationError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'IterationError: {0}'.format(self.message)

