
from lifty.sage_types import *
import random
import itertools as it
import sys

from lifty.kernel.triangulation import Segment, Edge, MarkedEdge, Vertex, Triangulation, LiftedTriangulation
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

        self.degree = len(self._coefficients)-1
        self._roots = None
        self._critical_pts = None
        self._postcritical_set = None
        self._pc_triangulation = None
        self._lifted_pcs = None
        self._lifted_tri = None
        self._derivative = None

    def coefficients(self):
        return self._coefficients
    
    def coeff(self, k):
        return self._coefficients[self.degree-k]

    def evaluate(self,z): # z should be an algebraic integer, if we want to return an algebraic integer
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
            self._critical_pts = [Qroots[i][0] for i in range(len(Qroots))]
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
            C = self.critical_pts()
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
                    return 'Error: either this polynomial is not post-critally finite, or max_iterations needs to be increased.'
                else:
                    continue
            self._postcritical_set = list(set([pcp for c in C for pcp in PCS[c]]))
        return self._postcritical_set

    def lifts(self,z):
        """Return all lifts of the point z, which must be an algebraic number.
        """
        c = self.coeff(0)-z
        coeffs = self.coefficients()[:-1]+[c]
        Q = TopPoly(coeffs)
        Qroots = Q.roots()
        lifts = set([Qroots[i][0] for i in range(len(Qroots))])
        return lifts

    def lifted_pcs(self):
        """Return lifts of the postcritical set. Results are cached.
        """
        if self._lifted_pcs == None:
            lifted_pcs = []
            for z in self.postcritical_set():
                lifts = self.lifts(z)
                for l in lifts:
                    lifted_pcs.append(l)
            self._lifted_pcs = lifted_pcs
        return self._lifted_pcs

    def dist_to_pcs(self,segment):
        """Returns minimal distance between the given edge and points of the 
            postcritical set which are not vertices of the edge. This should probably
            be moved out of the TopPoly class eventually.
        """
        a,b = CIF(segment[0]), CIF(segment[1])
        pcs = [p for p in self.postcritical_set() if p not in segment]
        m = Infinity
        for i in range(len(pcs)):
            p = CIF(pcs[i])
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
## max(|p'(z)|) <= max( |a_n*z^n| + |a_{n-1}*z^{n-1}| + ... + |a_1*z| + |a_0| ). If we maximize the latter
## quantity over the disk whose diameter is the segment A-->B, then the inequality still holds, and the 
## max is attained on the boundary of the disk. Such points are of the form z_0 + r_0*exp(i*theta), where 
## z_0 is the center of the disk and r_0 is the radius. Since |exp(i*theta)|=1, we can apply the binomial
## theorem to expand each term, then use the triangle inequality again to something that does not depend
## on theta, and his therefore constant. This is what the below function computes.

## Note that we do all of the above to avoid having to maximize a polynomial, at the cost of most likely 
## getting a worse bound on | ((p(B) - p(A))/(B - A) | than could have been obtained otherwise. My guess 
## is that this is worth it, but maybe there are improvements to be made here.



    def max_deriv(self,segment):
        """If Q(x) = a_n*z^n+a_{n-1}*z^{n-1}+ ... + a_1*z + a_0 is the derivative of self,
            then return the maximum of |a_n*z^n|+|a_{n-1}*z^{n-1}|+ ... + |a_1*z| + |a_0|
            on the disk whose diameter is segment. We really only care about the max of |Q(x)| on 
            segment, but it is easier computationally to take the max over a disk, since
            by the maximum modulus thm we know max(|Q(x)|) is attained on the boundary of the disk.
            Since max(|Q(x)|) <=  max(|a_n*z^n|+|a_{n-1}*z^{n-1}|+ ... + |a_1*z| + |a_0|), and 
            the latter is constant on the boundary of the disk, we simplify the computation
            by just returning that constant value, instead of actually maximizing |Q(x)|. 
            This should probably be moved out of the TopPoly class eventually.
        """
        Q = self.derivative()
        s1, s2 = CIF(segment[0]), CIF(segment[1])
        z0 = s1/2 + s2/2
        r0 = abs(z0-s1)
        d = Q.degree
        max_deriv = 0
        for i in range(d+1):
            a = CIF(Q.coeff(i))
            for k in range(d+1):
                max_deriv += abs(a)*binomial(i,k)*abs((z0**k)*(r0**(i-k)))
        return max_deriv


    def lift_edge(self,edge):
        """Return all lifts of edge, as a list of lists of algebraic numbers.
        """
        Q = self.derivative()
        a,b = edge.point(0), edge.point(1)
        a_lifts, b_lifts = self.lifts(a), self.lifts(b)
        
        lifted_edges = []

        ## to determine how many edges should have their initial point at each of the lifts of a,
        ## we lift a point c along the line from a to b that is very close to a, so that each lift of 
        ## c is distance <D from some lift of a. For a given lift A of a, the number of edges emanating 
        ## from A is the number of lifts of c that are D-close to A.
        t = QQ(1)/100000
        def init_edge(t, lifted_edges):
            c = b*t + (1-t)*a
            c_lifts = self.lifts(c)
            for A in a_lifts:
                for C in c_lifts:
                    D = self.dist_to_pcs([a,c])/(2*self.max_deriv([A,C]))
                    if abs(CIF(A)-CIF(C)) < D:
                        lifted_edges.append([A])
            return t, lifted_edges

        # if t is too large, we may not get all edges, so we keep making t smaller until we get self.degree edges
        while len(lifted_edges) < self.degree:
            t, lifted_edges = init_edge(t/10, lifted_edges)

        assert len(lifted_edges) == self.degree # this should always be true, but could fail if there is a bug or bad math.

        ## we extend the lift of each edge one segment at a time, until we have the full lift. Each
        ## iteration of below function extends by one segment. 
        def extend_edges(t,dt,interval,lifted_edges,finished):
            s1 = b*t+(1-t)*a
            s2 = b*(t+dt)+(1-(t+dt))*a
            lifted_segs = []
    
            seg = [s1,s2]
            s1_lifts = self.lifts(s1)
            s2_lifts = self.lifts(s2)
    
            for S1,S2 in it.product(s1_lifts,s2_lifts):
                
                ## if A and B are lifts of a and b, and abs(A-B) < D, then the image under self of the segment A--->B 
                ## is isotopic to the segment a--->b. So the segment A--->B is isotopic to a lift of a--->b.
                D = self.dist_to_pcs(seg)/(2*self.max_deriv([S1,S2]))
                if abs(CIF(S1)-CIF(S2)) < D:
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
                                if abs(CIF(lifted[-1]) - CIF(B)) < D:
                                    lifted.append(B)
                                    finished[i] = True
                                    break
                        ## if we are unable to jump to the end, find the lifted segment that should be joined to 
                        ## the edge in progress to extend it.            
                        if not finished[i]:
                            for i in range(len(lifted_segs)):
                                seg = lifted_segs[i]
                                if CIF(lifted[-1]).overlaps(CIF(seg[0])):
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
            for i in range(T.num_edges()):
                e = T.edge(i)
                e_lifted_edges = self.lift_edge(e)
                print('edge {} lifted.'.format(e.index()))

                for points in e_lifted_edges:
                    edge = MarkedEdge(points, i, None)
                    all_edges.append(edge)

            # set indices so they agree with the order in which edges are listed
            for i in range(len(all_edges)):
                e = all_edges[i]
                e._index = i

            ## check that any two distinct points in the lifted postcritical set are at least as far apart 
            ## as the diameter of their intervals. This ensures that we can test equality of points in lifted_pcs
            ## by comparing their intervals. This is important because confirming that two algebraic numbers 
            ## are equal is computationally hard, but checking that two intervals overlap is easy.
            max_pt = max([CIF(p).abs() for p in self.lifted_pcs()])
            max_dist = max([abs(CIF(a)-CIF(b)) for a,b in it.product(self.postcritical_set(),self.postcritical_set())])
            D = max_dist/(2*self.max_deriv([max_pt,-max_pt]))
            max_interval_diam = max([CIF(a).diameter() for a in self.lifted_pcs()])
            assert D > 2*max_interval_diam

            T = LiftedTriangulation(all_edges, self)
            self._lifted_tri = T
        return self._lifted_tri


#   def plot_lifted_edge(self, edge):
#       colors = ['red','blue','green']
#       lifted_edge = self.lift_edge(edge)
#       G = Graphics()
#       j=0
#       for e in lifted_edge:
#           G += line([CIF(e[i]) for i in range(len(e))],color=colors[j%3])
#           j+= 1
#       return G.show(dpi=1000)


    def pc_triangulation(self):
        """Return a triangulation whose vertices are the postcritical set. This will always be a 
            linear triangulation (each edge is a single line segment). Note that at this point 
            we don't include edges with endpoints at infinity. Results are cached.
        """
        if self._pc_triangulation == None:
            PCL = self.postcritical_set()
            L = len(PCL)
            edge_list = []
            for i in range(L-1):
                for j in range(i+1,L):
                    a,b = PCL[i],PCL[j]
                    e = Edge([a,b])
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
                                    if s.on_seg(PCL[k]):
                                        to_trash = True
                                        break
                        if not to_trash:
                            edge_list.append(e)

            for i in range(len(edge_list)):
                edge_list[i]._index = i

            T = Triangulation(edge_list, self)
            self._pc_triangulation = T

        return self._pc_triangulation


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

