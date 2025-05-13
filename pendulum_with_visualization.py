import networkx as nx
import numpy as np
from constrainthg.hypergraph import Hypergraph, Node
import constrainthg.relations as R
import matplotlib.pyplot as plt

#########################
###  Set Nodes ##########
#########################

hg = Hypergraph()

# -- Nodes
g = Node("gravity", -10)
r = Node("radius", 0.5)
theta0 = Node("theta0", 3.14159/6)
theta = Node("theta")
d_theta = Node("delta theta")
s_theta = Node("sine theta")
F = Node("gravitational force")
omega0 = Node("omega0", 0.0)
omega = Node("omega")
d_omega = Node("delta omega")
c = Node("damping coeff", 1.5)
alpha = Node("alpha")
time_step = Node("time_step", .03)
time = Node("time", 0.0)

########################################
##### Define Custom Relations ##########
########################################

# For complex systems, it is recommended to define these elsewhere and import
# See microgrid

def integrate(s1, s2, s3, **kwargs):
    return s1 * s3 + s2

def closed_form_theta(s1, s2, s3, s4, **kwargs):
      θ0 = kwargs['s1']
      g_ = kwargs['s2']
      ℓ  = kwargs['s3']
      t  = kwargs['s4']
      return θ0 * np.cos(np.sqrt(g_/ℓ) * t)


####################################
###### Define Edges ################
####################################

# Edges from Previous Pendulum Example, but they are updated to current
#best practice


hg.add_edge(theta0, 
            target= theta, 
            rel=R.Rmean, 
            label='theta0->theta')

hg.add_edge(omega0, 
            target=omega, 
            rel=R.Rmean, 
            label='omega0->omega')
# Note that the target (g/r) is not created in the initial node list.
# Targets and sources that are not located within the list are added. 

hg.add_edge({'s1': g, 's2': r}, 
            target='g/r', 
            rel=R.Rdivide, 
            label='(g,r)->b1')

# In case you are wondering why a dictionary is used here instead of a list or set,
# it is because the ordering of the arguments matters for this input equation
# g/r, and r/g are very different quantities. 


hg.add_edge(theta, 
            target=s_theta, 
            rel=R.Rsin, 
            label='theta->sine')
hg.add_edge({s_theta, 'g/r'}, 
            target=F, 
            rel=R.Rmultiply, 
            label='(sine, b1)->F')
hg.add_edge({omega, c}, 
            target='beta2', 
            rel=R.Rmultiply, 
            label='(omega, c)->b2')

# hg.addEdge(F, alpha, R.Rmean, label='F->alpha', index_offset=1)

hg.add_edge({'s1':F, 's2':'beta2'}, 
            target= alpha, 
            rel=R.Rsubtract, 
            label='(F, b2)->alpha', edge_props='LEVEL', index_offset=1)


# These following edges has been updated to utilize the updated, index_via function
hg.add_edge({'s1': alpha,'s2': omega,'s3': time_step}, 
            target=omega,
            rel=integrate, 
            label='int_w',
            index_via= lambda s1, s2, s3, **kwargs: s1 - 1 == s2)
            

hg.add_edge({'s1': omega,'s2': theta,'s3': time_step}, 
            target=theta, 
            rel=integrate, 
            label='(omega, theta, t)->theta',
            index_via=lambda s1, s2, s3, **kwargs: s1 - 1 == s2)

hg.add_edge({'s1':theta,'s2':omega}, 
            target='final theta', 
            rel=R.equal('s1'), 
            via=lambda s1, s2, **kwargs : abs(s1) < .05 and abs(s2) < .05, edge_props='LEVEL')


##############################################
###### Introduction of Global Time Tracking ##
##############################################


hg.add_edge({time, time_step},
    target=time,         
    rel=R.Rsum,       
    label='time',  
    index_offset=1)

hg.add_edge(
    {'s1': theta,
     's2': omega,
     's3': time},               
    target='final theta',
    rel=R.equal('s1'),
    via=lambda s1, s2, s3, **k: abs(s1)<.05 
                               and abs(s2)<.05
                               and s3 >= 0.0,
    edge_props='LEVEL')
# -------------------- 2) Symbol map --------------------

symbol_map = {
    'Rnull':         '0',
    'Rsum':          '+',
    'Rmultiply':     '×',
    'Rsubtract':     '–',
    'Rdivide':       '÷',
    'Rceiling':      '⌈⌉',
    'Rfloor':        '⌊⌋',
    'Rfloor_divide': '//',
    'Rnegate':       '−',
    'Rinvert':       '¹⁄x',
    'Rmean':         'μ',
    'Rmax':          'max',
    'Rmin':          'min',
    'Rsame':         '=',
    'Rall':          '∧',
    'Rany':          '∨',
    'Rxor':          '⊕',
    'Rnot_any':      '⊬',
    'Rnot':          '¬',
    'Rincrement':    '++',
    'Rfirst':        '1st',
    'Requal':        '=',    # equal()
    'Rcyclecounter': '≥',
    'Rsin':          'sin',
    'Rcos':          'cos',
    'Rtan':          'tan',
    'Rlist':         '[]',
    'Rtuple':        '()',
    'Rdict':         '{}',
    'integrate':     '∫',
}

# -------------------- 3) Solve & collect time series --------------------

t = hg.solve('final theta', to_print=False)
print(t)
# print(t.printTree())

getTimes = lambda l : [time_step.static_value * i for i in range(len(l))]
thetas = t.values[theta.label]
omegas = t.values[omega.label]
plt.plot(getTimes(thetas), thetas)
plt.plot(getTimes(omegas), omegas)
plt.legend(['theta', 'omega'])
plt.xlabel('Time (s)')
plt.ylabel('Rad, Rad/s')
plt.title('Pendulum Simulation')
plt.show()


def draw_hypergraph_bipartite(hg, symbol_map,
                              figsize=(7,6),
                              layout_fn=nx.spring_layout,
                              layout_kwargs=None):
    """
    Draw a CHG Hypergraph hg as a bipartite directed graph.
    
    Parameters
    ----------
    hg : Hypergraph
      Your ConstraintHG hypergraph.
    symbol_map : dict
      Maps relation-function names (edge.relation.__name__) to symbols.
      Falls back to edge.label if missing.
    figsize : tuple
      Matplotlib figure size.
    layout_fn : callable
      A NetworkX layout function (e.g. nx.spring_layout).
    layout_kwargs : dict
      kwargs to pass to layout_fn (e.g. {'k':1.0,'iterations':100}).
    
    Returns
    -------
    fig, ax : the matplotlib figure and axes.
    """
    # 1) Build the bipartite DiGraph
    B = nx.DiGraph()
    var_nodes = list(hg.nodes.keys())
    B.add_nodes_from(var_nodes, bipartite=0)
    
    for eid, edge in hg.edges.items():
        fnode = f"f:{eid}"
        # lookup relation name
        func = getattr(edge, 'relation', None) or getattr(edge, 'func', None)
        name = func.__name__ if func else None
        sym  = symbol_map.get(name, edge.label or '?')
        
        B.add_node(fnode, bipartite=1, label=sym)
        # sources → function node
        for src in edge.source_nodes.values():
            if isinstance(src, Node):
                B.add_edge(src.label, fnode)
        # function node → target
        if isinstance(edge.target, Node):
            B.add_edge(fnode, edge.target.label)
    
    # 2) Compute layout
    if layout_kwargs is None:
        layout_kwargs = {'k':1.0, 'iterations':100, 'seed':42}
    pos = layout_fn(B, **layout_kwargs)
    
    # 3) Prepare node labels
    node_labels = {
        n: (n if B.nodes[n]['bipartite']==0 else B.nodes[n]['label'])
        for n in B.nodes
    }
    
    # 4) Draw
    fig, ax = plt.subplots(figsize=figsize)
    
    # 4a) nodes (squares) at zorder=1
    nodes = nx.draw_networkx_nodes(
        B, pos,
        node_shape='s',
        node_size=600,
        node_color=[
            'skyblue' if B.nodes[n]['bipartite']==0 else 'gray'
            for n in B.nodes
        ],
        ax=ax
    )
    nodes.set_zorder(1)
    
    # 4b) edges (arrows) at zorder=2
    edges = nx.draw_networkx_edges(
        B, pos,
        arrows=True,
        arrowstyle='-|>',
        arrowsize=12,
        width=0.5,
        edge_color='gray',
        ax=ax
    )
    for art in edges:
        art.set_zorder(2)
    
    # 4c) labels at zorder=3
    texts = nx.draw_networkx_labels(
        B, pos,
        labels=node_labels,
        font_size=10,
        ax=ax
    )
    for txt in texts.values():
        txt.set_zorder(3)
    
    ax.set_title("Hypergraph as Bipartite Directed Graph\n(vars in blue, funcs in gray)")
    ax.axis('off')
    plt.tight_layout()
    return fig, ax


fig, ax = draw_hypergraph_bipartite(hg, symbol_map)
plt.show()