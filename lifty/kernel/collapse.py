import copy





def collapse_edge(CTri, edge_index):
    '''collapse the edge e with index edge_index. The square spanned by e consists of the two
        triangles adjacent to it, and each of these triangles is collapsed by identifying 
        its other two edges, and the two vertices that are the endpoint of e. 

        TO DO: ArcSystems associated to this triangulation are adjusted as needed.

        Warning: this operation changes the triangulation.
    '''

    T = CTri
    assert T.is_collapsible(edge_index), 'this edge is not collapsible'

    # place the edge so that it's vertical with the postcritical vertex at bottom.

    e = T.edge(edge_index)
    

    if e.vertex(0).is_ideal():
        v_bot = e.vertex(0)
        v_top = e.vertex(1)
        spokes = v_top.incident_edges()
        e_index = spokes.index(e.negative_rep())

        # get the (index of the) signed bottom-left and bottom-right edges of the square, with 
        # the sign relative to incidence at v_bot
        e_bot_left_index = (v_bot.incident_edges().index(e.positive_rep())+1)%v_bot.degree()
        e_bot_right_index = (v_bot.incident_edges().index(e.positive_rep())-1)%v_bot.degree()

    elif e.vertex(1).is_ideal():
        v_bot = e.vertex(1)
        v_top = e.vertex(0)
        spokes = v_top.incident_edges()
        e_index = spokes.index(e.positive_rep())
            
        e_bot_left_index = (v_bot.incident_edges().index(e.negative_rep())+1)%v_bot.degree()
        e_bot_right_index = (v_bot.incident_edges().index(e.negative_rep())-1)%v_bot.degree()

    # signed bottom-left and bottom-right edges
    e_bot_left = v_bot.incident_edges()[e_bot_left_index]
    e_bot_right = v_bot.incident_edges()[e_bot_right_index]


    degree = v_top.degree()

    # cyclically re-order the spoke edges so that e.signed_rep(sign) is the 0th edge, where
    # sign=1 if e is oriented downward, and sign=-1 otherwise.
    spokes = [spokes[(e_index+i)%degree] for i in range(len(spokes))]

    # get the triangles that meet v_top, ordered clockwise starting with the triangle to the left of e.
    incident_triangles = [spoke.adjacent_triangle() for spoke in spokes]


    # remove the two triangles adjacent to e from the triangle collection for this triangulation.
    triangles = list(T._triangles)
    triangles.remove(incident_triangles[0])
    triangles.remove(incident_triangles[-1])

    # remove the edge e, and the upper-left and upper-right edge of the square spanned by e.
    edges = list(T._edges)
    edges.remove(spokes[0].parent_edge())
    if degree > 2:
        edges.remove(spokes[1].parent_edge())
        edges.remove(spokes[-1].parent_edge())
    elif degree == 2:
        edges.remove(e_bot_left.parent_edge())
        edges.remove(e_bot_right.parent_edge())

    # remove the vertex v_top
    vertices = list(T._vertices)
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
    if degree == 2:
        edge_vertices = list(spokes[1].parent_edge()._vertices)
        if spokes[1].sign() == 1:
            edge_vertices[0] = v_bot
        else:
            edge_vertices[1] = v_bot
        spokes[1].parent_edge()._vertices = tuple(edge_vertices)


    if degree > 2:
        # now replace triangle edges for the upper-right triangle, which is first after the removed
        # triangles in the cyclic ordering around v_top 
        t_first = incident_triangles[1]
        tri_signed_edges = list(t_first._signed_edges)
        index = tri_signed_edges.index(spokes[1])
        tri_signed_edges[index] = e_bot_right

        # finalize these edge replacements
        t_first._signed_edges = tuple(tri_signed_edges)

        # now replace triangle edges for the upper-left triangle, which is last before the removed
        # triangles in the cyclic ordering around v_top 
        t_last = incident_triangles[-2]
        tri_signed_edges = list(t_last._signed_edges)
        index = tri_signed_edges.index(spokes[-1].opp_signed_edge())
        tri_signed_edges[index] = e_bot_left.opp_signed_edge()

        # finalize these edge replacements
        t_last._signed_edges = tuple(tri_signed_edges)


        # fix adjacent triangles for signed edges
        e_bot_left.opp_signed_edge()._adjacent_triangle = incident_triangles[-2]
        e_bot_right._adjacent_triangle = incident_triangles[1]

    elif degree == 2:

        t_right = e_bot_right.opp_signed_edge().adjacent_triangle()
        tri_signed_edges = list(t_right._signed_edges)
        index = tri_signed_edges.index(e_bot_right.opp_signed_edge())
        tri_signed_edges[index] = spokes[1].opp_signed_edge()

        t_right._signed_edges = tuple(tri_signed_edges)

        t_left = e_bot_left.adjacent_triangle()
        tri_signed_edges = list(t_left._signed_edges)
        index = tri_signed_edges.index(e_bot_left)
        tri_signed_edges[index] = spokes[1]

        t_left._signed_edges = tuple(tri_signed_edges)


        # fix adjacent triangles for signed edges
        spokes[1]._adjacent_triangle = t_left
        spokes[1].opp_signed_edge()._adjacent_triangle = t_right

    # now fix indices and finalize changes for triangulation
    for i in range(len(vertices)):
        v = vertices[i]
        v._index = i
        v._degree = None # this will be recomputed and cached next time v.degree() is called.
        v._incident_edges = None # this will be recomputed and cached next time v.incident_edges() is called.
    T._vertices = tuple(vertices)

    for i in range(len(edges)):
        e = edges[i]
        e._index = i
        e._positive_rep._index = i
        e._negative_rep._index = i
    T._edges = tuple(edges)

    for i in range(len(triangles)):
        t = triangles[i]
        t._index = i
    T._triangles = tuple(triangles)

    # make sure the resulting triangulation makes sense
    labels = [triangle.signed_edge(i).label() for triangle in T.triangles() for i in range(3)]
    for label in labels:
        assert labels.count(label) == 1
        assert ~int(label) in labels

    return True

def collapse_to_ideal(CTri):
    T = CTri
    T_copy = copy.deepcopy(T)

    collapsible = T_copy.collapsible_edges()

    while len(collapsible) > 0:
        _ = T_copy.collapse_edge(collapsible[0].index())
        collapsible = T_copy.collapsible_edges()

    assert all([v.is_ideal() for v in T_copy.vertices()])

    return T_copy












