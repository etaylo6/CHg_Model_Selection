from constrainthg.hypergraph import Hypergraph, Node
import constrainthg.relations as R
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

ss=.03
pivot1_x, pivot1_y = -1.0, 0.5  
pivot2_x, pivot2_y =  1.2, 0.8  
#########################
### Nodes ##########
#########################

hg = Hypergraph()

# -- Nodes
g_1 = Node("gravity_1", -10)
r_1 = Node("radius_1", 0.5)
theta0_1 = Node("theta0_1", 3.14159/6)
theta_1 = Node("theta_1")
d_theta_1 = Node("delta theta_1")
s_theta_1 = Node("sine theta_1")
F_1 = Node("gravitational force_1")
omega0_1 = Node("omega0_1", 0.0)
omega_1 = Node("omega_1")
d_omega_1 = Node("delta omega_1")
c_1 = Node("damping coeff_1", 1.5)
alpha_1 = Node("alpha_1")
time_step_1 = Node("time_step_1", ss)


##########################
##### Relations ##########
##########################

# For complex systems, it is recommended to define these elsewhere and import
# See microgrid

def integrate(s1, s2, s3, **kwargs):
    return s1 * s3 + s2



#PENDULUM 1#


############################
######Edges ################
############################

# Edges from Previous Pendulum Example, but they are updated to current
#best practice


hg.add_edge(theta0_1, 
            target= theta_1, 
            rel=R.Rmean, 
            label='theta0->theta')

hg.add_edge(omega0_1, 
            target=omega_1, 
            rel=R.Rmean, 
            label='omega0->omega')
# Note that the target (g/r) is not created in the initial node list.
# Targets and sources that are not located within the list are added. 

hg.add_edge({'s1': g_1, 's2': r_1}, 
            target='g/r_1', 
            rel=R.Rdivide, 
            label='(g,r)->b1')

# In case you are wondering why a dictionary is used here instead of a list or set,
# it is because the ordering of the arguments matters for this input equation
# g/r, and r/g are very different quantities. 


hg.add_edge(theta_1, 
            target=s_theta_1, 
            rel=R.Rsin, 
            label='theta->sine')
hg.add_edge({s_theta_1, 'g/r_1'}, 
            target=F_1, 
            rel=R.Rmultiply, 
            label='(sine, b1)->F')
hg.add_edge({omega_1, c_1}, 
            target='beta2_1', 
            rel=R.Rmultiply, 
            label='(omega, c)->b2')

# hg.addEdge(F, alpha, R.Rmean, label='F->alpha', index_offset=1)

hg.add_edge({'s1':F_1, 's2':'beta2_1'}, 
            target= alpha_1, 
            rel=R.Rsubtract, 
            label='(F, b2)->alpha', edge_props='LEVEL', index_offset=1)


# These following edges has been updated to utilize the updated, index_via function
hg.add_edge({'s1': alpha_1,'s2': omega_1,'s3': time_step_1}, 
            target=omega_1,
            rel=integrate, 
            label='(alpha, omega, t)->omega',
            index_via= lambda s1, s2, s3, **kwargs: s1 - 1 == s2)
            

hg.add_edge({'s1': omega_1,'s2': theta_1,'s3': time_step_1}, 
            target=theta_1, 
            rel=integrate, 
            label='(omega, theta, t)->theta',
            index_via=lambda s1, s2, s3, **kwargs: s1 - 1 == s2)

# The path to 'final theta' is only viable if the position and velocity of the pendulum are less than a threshold
hg.add_edge({'s1':theta_1,'s2':omega_1}, 
            target='final theta_1', 
            rel=R.equal('s1'), 
            via=lambda s1, s2, **kwargs : abs(s1) < .05 and abs(s2) < .05, edge_props='LEVEL')

###################### Pendulum #2#####################################


# -- Nodes
g = Node("gravity", -10)
r = Node("radius", 1)
theta0 = Node("theta0", 3.14159/2)
theta = Node("theta")
d_theta = Node("delta theta")
s_theta = Node("sine theta")
F = Node("gravitational force")
omega0 = Node("omega0", 0.0)
omega = Node("omega")
d_omega = Node("delta omega")
c = Node("damping coeff", 1.5)
alpha = Node("alpha")
time_step = Node("time_step", ss)

time = Node("t")
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
            label='(alpha, omega, t)->omega',
            index_via= lambda s1, s2, s3, **kwargs: s1 - 1 == s2)
            

hg.add_edge({'s1': omega,'s2': theta,'s3': time_step}, 
            target=theta, 
            rel=integrate, 
            label='(omega, theta, t)->theta',
            index_via=lambda s1, s2, s3, **kwargs: s1 - 1 == s2)


# The path to 'final theta' is only viable if the position and velocity of the pendulum are less than a threshold
hg.add_edge({'s1':theta,'s2':omega}, 
            target='final theta', 
            rel=R.equal('s1'), 
            via=lambda s1, s2, **kwargs : abs(s1) < .05 and abs(s2) < .05, edge_props='LEVEL')



##################################
#### Combined Viability Condition#
##################################

hg.add_edge({'s1' : 'final theta_1', 's2' : 'final theta'},
            target='check',
            rel=R.equal('s1'), 
            via=lambda s1, s2, **kwargs : abs(s1) < .05 and abs(s2) < .05)



t = hg.solve('check', to_print=False)
print(t)
# print(t.printTree())

thetas1 = t.values[theta_1.label]
thetas2 = t.values[theta.label]
times   = np.arange(len(thetas1)) * ss

# Compute (x,y) for each pendulum **including** pivot offsets
x1 = pivot1_x + r_1.static_value * np.sin(thetas1)
y1 = pivot1_y - r_1.static_value * np.cos(thetas1)
x2 = pivot2_x + r.static_value   * np.sin(thetas2)
y2 = pivot2_y - r.static_value   * np.cos(thetas2)

# Set up plot
fig, ax = plt.subplots()
ax.set_aspect('equal', 'box')

# Compute bounds to include pivots + swings
all_x = np.hstack(([pivot1_x, pivot2_x], x1, x2))
all_y = np.hstack(([pivot1_y, pivot2_y], y1, y2))
margin = max(r_1.static_value, r.static_value) * 0.2

ax.set_xlim(all_x.min() - margin, all_x.max() + margin)
ax.set_ylim(all_y.min() - margin, all_y.max() + margin)

ax.set_xlabel('x position (m)')
ax.set_ylabel('y position (m)')
ax.set_title('Two-Pendulum Animation')

# Prepare two rods and two bobs
rod1, = ax.plot([], [], lw=2)
bob1, = ax.plot([], [], 'o', markersize=12)
rod2, = ax.plot([], [], lw=2, linestyle='--')
bob2, = ax.plot([], [], 'o', markersize=12)
time_text = ax.text(0.05, 0.90, '', transform=ax.transAxes)

def init():
    for artist in (rod1, bob1, rod2, bob2, time_text):
        if hasattr(artist, 'set_data'):
            artist.set_data([], [])
    time_text.set_text('')
    return rod1, bob1, rod2, bob2, time_text

def update(i):
    # Pendulum 1
    rod1.set_data([pivot1_x, x1[i]], [pivot1_y, y1[i]])
    bob1.set_data(x1[i], y1[i])
    # Pendulum 2
    rod2.set_data([pivot2_x, x2[i]], [pivot2_y, y2[i]])
    bob2.set_data(x2[i], y2[i])
    time_text.set_text(f'Time = {times[i]:.2f} s')
    return rod1, bob1, rod2, bob2, time_text

ani = FuncAnimation(fig, update, frames=len(times),
                    init_func=init, blit=True,
                    interval=ss*1000)

plt.show()