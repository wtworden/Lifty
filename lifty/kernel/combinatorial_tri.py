

class CTriangulation():
    def __init__(self, triangles, vertices, vertex_markings):
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

        self._arc_systems = []

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

    def arc_systems(self):
        return self._arc_systems

    def arc_system(self,index):
        return self._arc_systems[index]

    def __repr__(self):
        return str(self._triangles)

    def collapsible_edges(self):
        return [e for e in self.edges() if self.is_collapsible(e.index())]

    def is_flippable(self, edge_index):
        e = self.edge(edge_index)
        if e.vertex(0).degree() > 1 and e.vertex(1).degree() > 1:
            return True

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

    def collapse_edge(self, edge_index):
        '''collapse the edge e with index edge_index. The square spanned by e consists of the two
            triangles adjacent to it, and each of these triangles is collapsed by identifying 
            its other two edges, and the two vertices that are the endpoint of e.

            ArcSystems associated to this triangulation are adjusted as needed.

            Warning: this operation changes the triangulation.
        '''
        print(self)
        print('edge:{}'.format(edge_index))
        assert self.is_collapsible(edge_index), 'this edge is not collapsible'

        # place the edge so that it's vertical with the postcritical vertex at bottom.
        e = self.edge(edge_index)
        

        if e.vertex(0).is_postcritical():
            v_bot = e.vertex(0)
            v_top = e.vertex(1)
            spokes = v_top.incident_edges()
            e_index = spokes.index(e.negative_rep())

            # get the (index of the) signed bottom-left and bottom-right edges of the square, with 
            # the sign relative to incidence at v_bot
            e_bot_left_index = (v_bot.incident_edges().index(e.positive_rep())-1)%v_bot.degree()
            e_bot_right_index = (v_bot.incident_edges().index(e.positive_rep())+1)%v_bot.degree()

        elif e.vertex(1).is_postcritical():
            v_bot = e.vertex(1)
            v_top = e.vertex(0)
            spokes = v_top.incident_edges()
            e_index = spokes.index(e.positive_rep())
        
            e_bot_left_index = (v_bot.incident_edges().index(e.negative_rep())-1)%v_bot.degree()
            e_bot_right_index = (v_bot.incident_edges().index(e.negative_rep())+1)%v_bot.degree()

        # signed bottom-left and bottom-right edges
        e_bot_left = v_bot.incident_edges()[e_bot_left_index]
        e_bot_right = v_bot.incident_edges()[e_bot_right_index]

        print('e_bot_left/right:{},{}'.format(e_bot_left,e_bot_right))

        degree = v_top.degree()

        # cyclically re-order the spoke edges so that e.signed_rep(sign) is the 0th edge, where
        # sign=1 if e is oriented downward, and sign=-1 otherwise.
        spokes = [spokes[(e_index+i)%degree] for i in range(len(spokes))]

        # get the triangles that meet v_top, ordered clockwise starting with the triangle to the left of e.
        incident_triangles = [spoke.adjacent_triangle() for spoke in spokes]


        # remove the two triangles adjacent to e from the triangle collection for this triangulation.
        triangles = list(self._triangles)
        triangles.remove(incident_triangles[0])
        triangles.remove(incident_triangles[-1])

        # remove the edge e, and the upper-left and upper-right edge of the square spanned by e.
        edges = list(self._edges)
        edges.remove(spokes[0].parent_edge())
        edges.remove(spokes[1].parent_edge())
        if degree > 2:
            edges.remove(spokes[-1].parent_edge())
        elif degree == 2:
            edges.remove(e_bot_right.parent_edge())

        # remove the vertex v_top
        vertices = list(self._vertices)
        vertices.remove(v_top)

        #now we need to clean things up, and fix vertices and edges of triangles.

        # first replace the vertex v_top of the incident triangles with v_bot
        for tri in incident_triangles[1:-1]: # if degree=2 then this list is empty, as it should be
            tri_vertices = list(tri._vertices)
            v_index = tri_vertices.index(v_top)
            tri_vertices[v_index] = v_bot
            tri._vertices = tuple(tri_vertices)

        # fix edge vertices for edges that had v_top as a vertex
        for signed_edge in spokes[2:-1]: # if degree<=3 then this list is empty, as it should be
            edge_vertices = list(signed_edge.parent_edge()._vertices)
            if signed_edge.sign() == 1:
                edge_vertices[0] = v_bot
            else:
                edge_vertices[1] = v_bot
            signed_edge.parent_edge()._vertices = tuple(edge_vertices)


        if degree > 2:
            # now replace triangle edges for the upper-left triangle, which is first after the removed
            # triangles in the cyclic ordering around v_top 
            t_first = incident_triangles[1]
            tri_signed_edges = list(t_first._signed_edges)
            index = tri_signed_edges.index(spokes[1])
            tri_signed_edges[index] = e_bot_left
    
            # finalize these edge replacements
            t_first._signed_edges = tuple(tri_signed_edges)
    
            # now replace triangle edges for the upper-right triangle, which is last before the removed
            # triangles in the cyclic ordering around v_top 
            t_last = incident_triangles[-2]
            tri_signed_edges = list(t_last._signed_edges)
            index = tri_signed_edges.index(spokes[-1].opp_signed_edge())
            tri_signed_edges[index] = e_bot_right.opp_signed_edge()
    
            # finalize these edge replacements
            t_last._signed_edges = tuple(tri_signed_edges)


            # fix adjacent triangles for signed edges
            e_bot_left._adjacent_triangle = incident_triangles[1]
            e_bot_right.opp_signed_edge()._adjacent_triangle = incident_triangles[-2]

        elif degree == 2:
            t_right = e_bot_right.adjacent_triangle()

            tri_signed_edges = list(t_right._signed_edges)
            index = tri_signed_edges.index(e_bot_right)
            tri_signed_edges[index] = e_bot_left

            t_right._signed_edges = tuple(tri_signed_edges)

            # fix adjacent triangles for signed edges
            e_bot_left._adjacent_triangle = t_right

        # now fix indices and finalize changes for triangulation
        for i in range(len(vertices)):
            v = vertices[i]
            v._index = i
            v._degree = None # this will be recomputed and cached next time v.degree() is called.
            v._incident_edges = None # this will be recomputed and cached next time v.incident_edges() is called.
        self._vertices = tuple(vertices)

        for i in range(len(edges)):
            e = edges[i]
            e._index = i
            e._positive_rep._index = i
            e._negative_rep._index = i
        self._edges = tuple(edges)

        for i in range(len(triangles)):
            t = triangles[i]
            t._index = i
        self._triangles = tuple(triangles)

        # make sure the resulting triangulation makes sense
        labels = [triangle.signed_edge(i).label() for triangle in self.triangles() for i in range(3)]
        for label in labels:
            assert labels.count(label) == 1
            assert ~int(label) in labels

    def collapse_to_ideal(self):
        collapsible = self.collapsible_edges()

        while len(collapsible) > 0:
            self.collapse_edge(collapsible[0].index())
            collapsible = self.collapsible_edges()



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
        return self._vertices[edge_index]


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

    def adjacent_tri_right(self):
        return self.positive_rep().adjacent_triangle()

    def adjacent_tri_left(self):
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

    def triangulation(self):
        return self._triangulation

    def index(self):
        return self._index

    def incident_edges(self):
        '''return edge incident to this vertex, ordered clockwise starting
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

    def is_postcritical(self):
        return self._is_postcritical

    def is_ideal(self):
        return self._is_ideal

    def degree(self):
        if self._degree == None:
            self._degree = len(self.incident_edges())
        return self._degree

    def __repr__(self):
        return str((self._incident_edges, self.is_postcritical()))

class Arc():
    def __init__(self, edge_intersections):
        self._vector = tuple(edge_intersections)

    def edge_weight(self,edge_index):
        return self._vector[edge_index]

    def vector(self):
        return self._vector

class ArcSystem():
    def __init__(self, arcs):
        self._arcs = arcs

    def arcs(self):
        return self._arcs

    def arc(self,index):
        return self._arcs[index]















