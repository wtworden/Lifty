import copy

def find_path(T1, T2):
    """Find a path of pachner flips that take the Curver triangulation T1 to the Curver Triangulation T2."""

    ### TO DO: everything

    pass



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

    ### TO DO: keep track of multi-arcs

    T = CTri

    assert T.is_flippable(edge_index)

    e = T.edge(edge_index)
    
    e_pos = e.positive_rep()
    e_neg = e.negative_rep()
    
    t1 = e_pos.adjacent_triangle()
    t2 = e_neg.adjacent_triangle()
    
    e_pos_index = t1.signed_edges().index(e_pos)
    e_neg_index = t2.signed_edges().index(e_neg)

    e1 = t1.signed_edge((e_pos_index+2)%3)
    e2 = t1.signed_edge((e_pos_index+1)%3)
    e3 = t2.signed_edge((e_neg_index+2)%3)
    e4 = t2.signed_edge((e_neg_index+1)%3)


    v1 = e.vertex(0)
    v3 = e.vertex(1)

    v2 = t1.opp_vert(e)
    v4 = t2.opp_vert(e)

    # flip e
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

    # make sure the resulting triangulation makes sense
    labels = [triangle.signed_edge(i).label() for triangle in T.triangles() for i in range(3)]
    for label in labels:
        assert labels.count(label) == 1
        assert ~int(label) in labels

    return True
