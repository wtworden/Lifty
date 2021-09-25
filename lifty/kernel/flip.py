import copy


    # the below drawing shows what happens when an oriented edge e is flipped. We also denote the new 
    # edge by e. The left/right triangles t1 and t2 are replaced by the top/bottom triangles
    # also denoted t1 and t2. The labeling of simplices here matches that used in flip_edge().
    # Note that it could happen that ei and ej are identified, for some pair ei,ej in [e1,e2,e3,e4],
    # but this will not have any affect on the edge flip


    #             v3                                 v3                                                 
    #             o                                  o                                                  
    #            /|\                                / \                                                     
    #       e2 /  |  \ e3                      e2 / t1  \ e3                                           
    #        /    |    \                        /         \                                                
    #      /      |e     \                    /    e        \                                                
    #  v2 o   t1  |   t2  o v4   ---->    v2 o------->-------o v4                                                          
    #      \      ^      /                    \             /                                                     
    #        \    |    /                        \    t2   /                                                               
    #       e1 \  |  / e4                      e1 \     / e4                                                       
    #            \|/                                \ /                                                         
    #             o                                  o                                                      
    #             v1                                 v1                                                     
    #

def flip_edge(CTri, edge_index):

    T = CTri

    # make sure this is a flippable edge
    assert T.is_flippable(edge_index), "edge is not flippable!"

    # the edge to flip
    e = T.edge(edge_index)
    
    # the two signed edges associated to e
    e_pos = e.positive_rep()
    e_neg = e.negative_rep()
    
    #the two triangles adjacent to e
    t1 = e_pos.adjacent_triangle()
    t2 = e_neg.adjacent_triangle()
    
    # the index (0,1, or 2) with respect to each triangle of the signed edge
    e_pos_index = t1.signed_edges().index(e_pos)
    e_neg_index = t2.signed_edges().index(e_neg)

    # the other two edges of t1 (see diagram above)
    e1 = t1.signed_edge((e_pos_index+2)%3)
    e2 = t1.signed_edge((e_pos_index+1)%3)
    # the other two edges of t2 (see diagram above)
    e3 = t2.signed_edge((e_neg_index+2)%3)
    e4 = t2.signed_edge((e_neg_index+1)%3)

    #the bottom and top vertices of the quad
    v1 = e.vertex(0)
    v3 = e.vertex(1)

    # the left and right vertices of the quad
    v2 = t1.opp_vert(e)
    v4 = t2.opp_vert(e)

    # adjust multi-arc intersection info
    for key in T.multi_arcs():
        M = T.multi_arc(key)
        for arc in M.arcs():
            # we will adjust the intersection_sequence for the arc, then recompute the geometric
            # and algebraic intersection vectors. 

            # first let's make sure this arc has no bigons:
            arc._remove_bigons()
            arc._remove_end_bigons()

            # first set the algebraic and geometric intersections to 0.
            alg = 0
            geom = 0

            # it will be useful to have a shorthand for the last index in the intersection_sequence
            last_index = len(arc.intersection_sequence())-1

            # now look through the intersection sequence for patterns that indicate that the new edge
            # will have an intersection. Basically, we just exhaust all possible sub-arcs types that
            # would intersect e post-flipping.
            for i in range(len(arc.intersection_sequence())):
                edge = arc.intersection_sequence()[i]
                if i == 0:
                    if arc.vertex(0) == v3:
                        if edge == e1 or edge == e4:
                            alg += 1
                            geom += 1
                    elif arc.vertex(0) == v1:
                        if edge == e2 or edge == e3:
                            alg += -1
                            geom += 1
                if edge == ~e1:
                    if (i == last_index) or (arc.intersection_sequence()[i+1] == e2) or (i+1 != last_index and arc.intersection_sequence()[i+2] == e3):
                        alg += -1
                        geom += 1
                if edge == ~e2:
                    if (i == last_index) or (arc.intersection_sequence()[i+1] == e1) or (i+1 != last_index and arc.intersection_sequence()[i+2] == e4):
                        alg += 1
                        geom += 1                    
                if edge == ~e3:
                    if (i == last_index) or (arc.intersection_sequence()[i+1] == e4) or (i+1 != last_index and arc.intersection_sequence()[i+2] == e1):
                        alg += 1
                        geom += 1
                if edge == ~e4:
                    if (i == last_index) or (arc.intersection_sequence()[i+1] == e3) or (i+1 != last_index and arc.intersection_sequence()[i+2] == e2):
                        alg += -1
                        geom += 1    
            arc._algebraic[e.index()] = alg
            arc._geometric[e.index()] = geom

    # now we flip e
    e._vertices = tuple([v2, v4])

    t1._signed_edges = tuple([e_pos, e3, e2])
    t2._signed_edges = tuple([e_neg, e1, e4])

    t1._signs = tuple([1, e3.sign(), e2.sign()])
    t2._signs = tuple([-1, e1.sign(), e4.sign()])

    t1._vertices = tuple([v3, v2, v4])
    t2._vertices = tuple([v1, v4, v2])

    # reset vertex incident edges and vertex degree to None. These will be recomputed and 
    # cached the next time the methods are called.
    for v in [v1, v2, v3, v4]:
        v._incident_edges = None
        v._degree = None

    # reset adjacent triangles for affected signed edges to None. These will be recomputed
    # and chached the next time they're called.
    e1._adjacent_triangle = None
    e2._adjacent_triangle = None
    e3._adjacent_triangle = None
    e4._adjacent_triangle = None
    e_pos._adjacent_triangle = None
    e_neg._adjacent_triangle = None


    # now re-compute the intersection_sequences for all arcs
    for key in T.multi_arcs():
        M = T.multi_arc(key)
        for arc in M.arcs():    
            arc._compute_intersection_sequence()

    # make sure the resulting triangulation makes sense
    labels = [triangle.signed_edge(i).label() for triangle in T.triangles() for i in range(3)]
    for label in labels:
        assert labels.count(label) == 1
        assert ~int(label) in labels

    return True

def flip_to_pc_tri(CTri):
    assert not CTri.has_material_vertices()
    T = copy.deepcopy(CTri)

    pc_multi_arc = T.multi_arc('pc_tri')

    # we'll keep track of the sequence of edge flips we do here
    flips = []

    # get arcs that intersect T
    arcs = [a for a in pc_multi_arc.arcs() if sum(a.geometric()) > 0]

    while len(arcs) > 0:
        arc = arcs[0]
        while len(arc.intersection_sequence()) > 0:
            to_flip = arc.intersection_sequence()[0].index()
            flips.append(to_flip)
            _ = flip_edge(T, to_flip)
        arcs = [a for a in pc_multi_arc.arcs() if sum(a.geometric()) > 0]

    return T, flips
















