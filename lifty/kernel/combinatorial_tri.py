
import lifty.kernel.collapse as collapse
import lifty.kernel.flip as flip

from lifty.sage_types import flatten


class CTriangulation():
    def __init__(self, triangles, vertices, vertex_markings, vertex_mappings, is_lifted_tri):
        self._is_lifted_tri = is_lifted_tri
        self._triangles = tuple(triangles)

        for i in range(len(self._triangles)):
            t = self._triangles[i]
            t._index = i
            t._triangulation = self

        # edges of the triangulation are oriented edge classes. Each contains two reps, one for each adjacent
        # triangle. One is the positive rep, and is associated to the adjacent triangle whose clockwise orientation
        # agrees with the orientation of the edge, and the other is the negative rep, associated to the adjacent
        # triangle whose clockwise orientation is opposite to the edge orientation.
        self._edges = tuple(sorted([triangle.edge(i) for triangle in triangles for i in range(3) if triangle.signed_edge(i).sign()==1],key=lambda x:x.index()))
        for e in self._edges:
            e._triangulation = self

        self._vertices = tuple(vertices)

        for i in range(len(self._vertices)):
            v = self._vertices[i]
            v._triangulation = self
            v._is_postcritical = vertex_markings[i][0]
            v._is_ideal = vertex_markings[i][1]
            if self._is_lifted_tri:
                v._maps_to = vertex_mappings[i]

        self._multi_arcs = {}

    def is_lifted_tri(self):
        return self._is_lifted_tri

    def vertices(self):
        return self._vertices

    def vertex(self, index):
        return self._vertices[index]

    def num_edges(self):
        return len(self._edges)

    def edges(self):
        return self._edges

    def edge(self,index):
        return self._edges[index]

    def triangles(self):
        return self._triangles

    def triangle(self, index):
        return self._triangles[index]

    def self_glued(self, triangle_index):
        return len(set([edge.index() for edge in self.triangle(triangle_index).edges()])) < 3

    def multi_arcs(self):
        return self._multi_arcs

    def multi_arc(self,key):
        return self._multi_arcs[key]

    def __repr__(self):
        return str(self._triangles)

    def is_flippable(self, edge_index):
        e = self.edge(edge_index)
        if e.vertex(0).degree() > 1 and e.vertex(1).degree() > 1:
            return True

    def flippable_edges(self):
        return [e for e in self.edges() if self.is_flippable(e.index())]

    def flip_edge(self, edge_index):        
        return flip.flip_edge(self, edge_index)


    def is_collapsible(self,edge_index):
        '''an edge is collapsible if neither adjacent triangle is self glued, and the edge
            has one vertex that is postcritical, and one that is not.
        '''
        e = self.edge(edge_index)
        t1, t2 = e.adjacent_triangles()
        # first check that neither adjacent triangle has two of its edges identified
        if not ( self.self_glued(t1.index()) or self.self_glued(t2.index()) ):
            is_ideal = [e.vertex(0).is_ideal(), e.vertex(1).is_ideal()]
            # next check that the edge has exactly one ideal vertex
            if True in is_ideal and False in is_ideal:
                return True
        return False


    def collapsible_edges(self):
        return [e for e in self.edges() if self.is_collapsible(e.index())]

    def collapse_edge(self, edge_index):
        '''collapse the edge e with index edge_index. The square spanned by e consists of the two
            triangles adjacent to it, and each of these triangles is collapsed by identifying 
            its other two edges, and the two vertices that are the endpoint of e.

            MultiArcs associated to this triangulation are adjusted along the way.

            Warning: this operation changes the triangulation.
        '''
        return collapse.collapse_edge(self,edge_index)

    def collapse_to_ideal(self):
        '''Successively collapse all edges of the triangulation that have one vertex that is ideal, 
            and one that is material. Returns a triangulation having only ideal vertices (i.e., every 
            vertex is either a postcritical point or infinity). This does not change the original 
            triangulation, and returns a new triangulation.
        '''

        return collapse.collapse_to_ideal(self)


class CTriangle():
    def __init__(self, edges, signs):
        self._signs = tuple(signs)
        self._triangulation = None

        self._index = None

        vertices = [None, None, None]
        for i in range(3):
            e = edges[i]
            if self._signs[i] == 1:
                vertices[(i+1)%3] = e.vertex(0)
                e._positive_rep = SignedEdge(e,e.index(),1)
                e._positive_rep._adjacent_triangle = self
            elif self._signs[i] == -1:
                vertices[(i+1)%3] = e.vertex(1)
                e._negative_rep = SignedEdge(e,e.index(),-1)
                e._negative_rep._adjacent_triangle = self                
        
        self._vertices = tuple(vertices)
        self._signed_edges = tuple([edges[i].signed_rep(self._signs[i]) for i in range(3)])

    def index(self):
        return self._index

    def signed_edges(self):
        return self._signed_edges

    def signed_edge(self, index):
        return self._signed_edges[index]

    def edges(self):
        return [e.parent_edge() for e in self.signed_edges()]

    def edge(self, index):
        return self.signed_edge(index).parent_edge()

    def vertices(self):
        return self._vertices

    def vertex(self, index):
        return self._vertices[index]

    def opp_vert(self, edge):
        edge_index = self.edges().index(edge)
        return self.vertex(edge_index)

    def opp_signed_edge(self, vert):
        vert_index = self.vertices().index(vert)
        return self.signed_edge(vert_index)

    def __repr__(self):
        return str(self._signed_edges)

class CEdge():
    def __init__(self, index, vertices):
        self._index = index
        self._vertices = tuple(vertices)
        self._triangulation = None

        self._positive_rep = None
        self._negative_rep = None

    def signed_rep(self, sign):
        if sign == 1:
            return self._positive_rep
        else:
            return self._negative_rep

    def positive_rep(self):
        return self._positive_rep

    def negative_rep(self):
        return self._negative_rep

    def index(self):
        return self._index

    def vertices(self):
        return self._vertices

    def vertex(self, index):
        return self._vertices[index]

    def adjacent_triangles(self):
        return [self.adjacent_tri_right(), self.adjacent_tri_left()]

    def adjacent_tri_left(self):
        return self.positive_rep().adjacent_triangle()

    def adjacent_tri_right(self):
        return self.negative_rep().adjacent_triangle()

    def __repr__(self):
        return str(self.index())

class SignedEdge():
    def __init__(self, parent, index, sign):
        self._index = index
        self._adjacent_triangle = None

        self._sign = sign

        self._parent_edge = parent

    def opp_signed_edge(self):
        if self._sign == 1:
            return self._parent_edge.negative_rep()
        else:
            return self._parent_edge.positive_rep()

    def label(self):
        return self._index if self._sign==1 else ~int(self._index)

    def parent_edge(self):
        return self._parent_edge

    def sign(self):
        return self._sign

    def index(self):
        return self._index

    def vertices(self):
        return self.parent_edge().vertices()

    def vertex(self, index):
        return self.parent_edge().vertex(index)

    def adjacent_triangle(self):
        return self._adjacent_triangle

    def __repr__(self):
        return '~'*int(self.sign()==-1) + str(self.index())

class CVertex():
    def __init__(self, index):
        self._index = index
        self._is_postcritical = None
        self._is_ideal = None
        self._triangulation = None
        self._degree = None
        self._incident_edges = None
        self._maps_to = None

    def triangulation(self):
        return self._triangulation

    def index(self):
        return self._index

    def incident_edges(self):
        '''return edge incident to this vertex, ordered counter-clockwise starting
            with the edge of smallest index. If the edge is oriented 
            into the vertex, we return its negative_rep. If the edge is
            oriented out from the vertex, return its positive_rep.
        '''
        if self._incident_edges == None:
            T = self.triangulation()
            triangles = T.triangles()
            incident_edges = []
            edge_pairs = []
            for t in triangles:
                for i in range(3):
                    if self == t.vertex(i):
                        edge_pairs.append((t.signed_edge((i-1)%3), t.signed_edge((i+1)%3).opp_signed_edge()))
            edge_pairs.sort(key = lambda x: x[0].index())
            p = edge_pairs.pop(0)
            incident_edges.append(p[0])
            while len(edge_pairs) > 0:
                for i in range(len(edge_pairs)):
                    if edge_pairs[i][0] == p[1]:
                        p = edge_pairs.pop(i)
                        incident_edges.append(p[0])
                        break
            self._incident_edges = incident_edges

        return self._incident_edges

    def next_incident_edge(self, edge):
        this = self.incident_edges().index(edge)
        return self.incident_edges()[(this+1)%self.degree()]


    def maps_to(self):
        return self._maps_to

    def is_postcritical(self):
        return self._is_postcritical

    def is_ideal(self):
        return self._is_ideal

    def degree(self):
        if self._degree == None:
            self._degree = len(self.incident_edges())
        return self._degree

    def __repr__(self):
        return str((self.incident_edges(), self.is_postcritical()))

class Arc():
    def __init__(self, Ctriangulation, vertices, algebraic_ints, geometric_ints, intersection_sequence=None):
        
        self._triangulation = Ctriangulation
        self._vertices = tuple([self._triangulation.vertex(vertices[0]), self._triangulation.vertex(vertices[1])])
        assert self._vertices[1].is_ideal()
        assert self._vertices[0].is_ideal()
        
        if intersection_sequence != None:
            self._intersection_sequence = []
            for i in intersection_sequence:
                i = int(i)
                if i<0:
                    self._intersection_sequence.append(self._triangulation.edge(~i).negative_rep())
                else:
                    self._intersection_sequence.append(self._triangulation.edge(i).positive_rep())

            self._geometric = None
            self._algebraic = None
            self._recompute_intersection_vectors()
        else:
            self._geometric = geometric_ints
            self._algebraic = algebraic_ints


            # it will be useful to also encode the arc as a sequence of edge intersections (i.e., start
            # at self._vertices[0] and walk along the arc, record edges crossed (with sign) as you go). This 
            # makes it especially easy to remove bigons.
            self._intersection_sequence = []
            if sum(self._geometric) > 0:
                T = self._triangulation
                v0 = self._vertices[0]
                v1 = self._vertices[1]
                # first we need to find the initial edge crossed:
                for t in T.triangles():
                    if v0 in t.vertices():
                        e_opp = t.opp_signed_edge(v0)
                        e_left = t.edge((t.signed_edges().index(e_opp)+1)%3)
                        e_right = t.edge((t.signed_edges().index(e_opp)+2)%3)
                        a = self._geometric[e_opp.index()]
                        b = self._geometric[e_left.index()]
                        c = self._geometric[e_right.index()]
                        if a > b + c:
                            print(a,b,c)
                            self._intersection_sequence.append(e_opp)
                            exit_edge = e_opp
                            enter_edge = e_opp.opp_signed_edge()
                            next_triangle = enter_edge.adjacent_triangle()
                            position = b
                            print(exit_edge, enter_edge, next_triangle,position)
                            break
                while exit_edge != None:
                    t = next_triangle
                    e_right = t.signed_edge((t.signed_edges().index(enter_edge)+1)%3)
                    e_left = t.signed_edge((t.signed_edges().index(enter_edge)+2)%3)
                    a = self._geometric[enter_edge.index()]
                    b = self._geometric[e_left.index()]
                    c = self._geometric[e_right.index()]
                    top_arcs = (b+c-a)/2
                    b_a = b - max([0,top_arcs])
                    if top_arcs >= 0 or position != b_a:
                        if b_a == position:
                            self._intersection_sequence.append(e_right)
                            exit_edge = e_right
                            enter_edge = e_right.opp_signed_edge()
                            next_triangle = enter_edge.adjacent_triangle()
                            position = top_arcs
                        elif b_a < position:
                            self._intersection_sequence.append(e_right)
                            exit_edge = e_right
                            enter_edge = e_right.opp_signed_edge()
                            next_triangle = enter_edge.adjacent_triangle()
                            if top_arcs >= 0:
                                position = top_arcs + (position - b_a)
                            else:
                                position = position - b_a - 1
                        elif b_a > position:
                            self._intersection_sequence.append(e_left)
                            exit_edge = e_left
                            enter_edge = e_left.opp_signed_edge()
                            next_triangle = enter_edge.adjacent_triangle()
                            position = position
                    else:
                        exit_edge = None

    def intersection_sequence(self):
        return self._intersection_sequence

    def vertex(self, index):
        return self._vertices[index]

    def vertices(self):
        return self._vertices

    def edge_weight(self,edge_index):
        return self._geometric[edge_index]

    def algebraic(self):
        return self._algebraic

    def geometric(self):
        return self._geometric

    def triangulation(self):
        return self._triangulation

    # bigons correspond to subsequences of self._intersection_sequence of the form x,~x.
    def _has_bigons(self):

        for i in range(len(self._intersection_sequence)-1):
            if self._intersection_sequence[i].label() == ~self._intersection_sequence[i+1].label():
                return True, self._intersection_sequence[i], (i,i+1)
        return False, None, (None,None)

    #remove all bigons recursively
    def _remove_bigons(self):
        bigons, _, (i,j) = self._has_bigons()
        while bigons:
            _ = self._intersection_sequence.pop(j)
            _ = self._intersection_sequence.pop(i)
            bigons, _, (i,j) = self._has_bigons()

    # A normal arc with endpoint at a vertex v of triangle t should have its first intersection
    # on the edge across from v. If the first intersection is with an edge that shares the vertex
    # v, then this gives a bigon (with one vertex v, hence an end bigon). 
    def _has_end_bigons(self):
        if len(self._intersection_sequence) > 0:
            first_intersection = self._intersection_sequence[0].index() 
            last_intersection = self._intersection_sequence[-1].index()
            if self.vertex(0) in self.triangulation().edge(first_intersection).vertices():
                return True, first_intersection, 0 
            if self.vertex(1) in self.triangulation().edge(last_intersection).vertices():
                return True, last_intersection, -1
        return False, None, None

    def _remove_end_bigons(self):
        end_bigons, edge, i = self._has_end_bigons()
        while end_bigons:
            self._intersection_sequence.pop(i)
            end_bigons, edge, i = self._has_end_bigons()


    def isotope_through_vertex(self, edge_index, vertex_index):
        """
        """
        T = self._triangulation
        assert T.vertex(vertex_index).is_ideal() == False, "you cannot isotope over an ideal vertex"
        

        if self.edge_weight(edge_index) == 0: #if the arc doesn't intersect this edge, do nothing.
            pass
        else:
            e_unsigned = T.edge(edge_index)
            v = T.vertex(vertex_index)
            degree = v.degree()
            if e_unsigned.vertex(0) == v:
                w = e_unsigned.vertex(1)
                e = e_unsigned.positive_rep()
            else:
                w = e_unsigned.vertex(0)
                e = e_unsigned.negative_rep()
            e_rel_v = v.incident_edges().index(e)
    
            v_incident_sans_e = v.incident_edges()[e_rel_v+1:]+v.incident_edges()[0:e_rel_v]
    
            for i in range(len(self.intersection_sequence())):
                edge = self.intersection_sequence()[i]
                if edge == e:
                    new_subseq = [a.opp_signed_edge() for a in v_incident_sans_e]
                    _ = self._intersection_sequence.pop(i)
                    self._intersection_sequence.insert(i,new_subseq)
                elif edge == e.opp_signed_edge():
                    new_subseq = list(reversed(v_incident_sans_e))
                    _ = self._intersection_sequence.pop(i)
                    self._intersection_sequence.insert(i,new_subseq)
                else:
                    pass
            self._intersection_sequence = flatten(self._intersection_sequence)
    
            self._remove_bigons()
            self._remove_end_bigons()


    def _recompute_intersection_vectors(self):
        geometric = []
        algebraic = []
        T = self.triangulation()
        for i in range(T.num_edges()):
            pos = self.intersection_sequence().count(T.edge(i).positive_rep())
            neg = self.intersection_sequence().count(T.edge(i).negative_rep())
            geometric.append(pos+neg)
            algebraic.append(pos-neg)
        self._geometric = geometric
        self._algebraic = algebraic
        

    def __repr__(self):
        return 'arc with vertices: [{}, {}], intersections: {}'.format(self.vertex(0), self.vertex(1), self.algebraic())


class MultiArc():
    def __init__(self, arcs):
        self._arcs = arcs

    def arcs(self):
        return self._arcs

    def arc(self,index):
        return self._arcs[index]

    def isotope_through_vertex(self, edge_index, vertex_index):
        for arc in self._arcs:
            arc.isotope_through_vertex(edge_index, vertex_index)

    def _recompute_intersection_vectors(self):
        for arc in self._arcs:
            arc._recompute_intersection_vectors()

    def __repr__(self):
        out = 'multi-arc with arcs:\n'
        for i in range(len(self.arcs())):
            out += '{}: {}\n'.format(i,self.arc(i).__repr__())

        return out















