from SystemModels import Beam_Lattice
import numpy as np
import matplotlib.pyplot as plt
from Forces import stationary_input_shaping_force
from SystemSolvers import Newmark, Static
import matplotlib.pyplot as plt
import matplotlib.animation as animation

system = Beam_Lattice()

system.add_beam_edge(
    number_of_elements=1, 
    E_modulus=2.1*10**11, 
    shear_modulus=7.9*10**10,
    primary_moment_of_area=2.157*10**-8,
    secondary_moment_of_area=1.113*10**-8,
    torsional_constant=3.7*10**-8,
    density=7850, 
    cross_sectional_area=1.737*10**-4, 
    coordinates=[[0, 0, 0], [0, 0, 1]],
    edge_polar_rotation=0
)

system.add_beam_edge(
    number_of_elements=1, 
    E_modulus=2.1*10**11, 
    shear_modulus=7.9*10**10,
    primary_moment_of_area=2.157*10**-8,
    secondary_moment_of_area=1.113*10**-8,
    torsional_constant=3.7*10**-8,
    density=7850, 
    cross_sectional_area=1.737*10**-4, 
    coordinates=[0, 0, 2],
    vertex_IDs=1,
    edge_polar_rotation=0
)

system.add_beam_edge(
    number_of_elements=1, 
    E_modulus=2.1*10**11, 
    shear_modulus=7.9*10**10,
    primary_moment_of_area=2.157*10**-8,
    secondary_moment_of_area=1.113*10**-8,
    torsional_constant=3.7*10**-8,
    density=7850, 
    cross_sectional_area=1.737*10**-4, 
    coordinates=[0, 0, 3],
    vertex_IDs=2,
    edge_polar_rotation=0
)

system.add_beam_edge(
    number_of_elements=1, 
    E_modulus=2.1*10**11, 
    shear_modulus=7.9*10**10,
    primary_moment_of_area=2.157*10**-8,
    secondary_moment_of_area=1.113*10**-8,
    torsional_constant=3.7*10**-8,
    density=7850, 
    cross_sectional_area=1.737*10**-4, 
    coordinates=[0, 0, 4],
    vertex_IDs=3,
    edge_polar_rotation=0
)

system.add_beam_edge(
    number_of_elements=1, 
    E_modulus=2.1*10**11, 
    shear_modulus=7.9*10**10,
    primary_moment_of_area=2.157*10**-8,
    secondary_moment_of_area=1.113*10**-8,
    torsional_constant=3.7*10**-8,
    density=7850, 
    cross_sectional_area=1.737*10**-4, 
    coordinates=[0, 0, 5],
    vertex_IDs=4,
    edge_polar_rotation=0
)

system.fix_vertices((0,))
system.add_forces({4: lambda t: [0, np.sin(t+2), 0, 0, 0, 0],
                   3: lambda t: [0, np.sin(t), 0, 0, 0, 0]})
np.set_printoptions(linewidth=1000, precision=2)

stationary_vertex_IDs = 1
input_vertex_IDs = 2
time_increment = 0.01
end_time = 30
scaling_factor = 500

stationary_input_shaping_force(system, np.arange(0, 2*np.pi, time_increment), stationary_vertex_IDs, input_vertex_IDs)

plt.scatter(np.arange(0, end_time, time_increment), [system.graph.vs[2]['force'](time)[3] for time in np.arange(0, end_time, time_increment)])

plt.show()

static_solver = Static(system).solve()
solver = Newmark(system, static_solver, end_time, time_increment).solve()
displaced_shape_points = solver.get_displaced_shape_position(scaling_factor)

plot_lines = list()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for displaced_shape_point in displaced_shape_points:
    plot_lines.append(ax.plot(displaced_shape_point[0, :, 0], displaced_shape_point[0, :, 1], displaced_shape_point[0, :, 2], linewidth=2.0))

def update(frame):
    
    for i, lines in enumerate(plot_lines):
        for line in lines:
            line.set_data_3d(*displaced_shape_points[i][frame, :].T)
    ax.set_title(f"t = {frame*time_increment:.3f} s")
    return plot_lines

ani = animation.FuncAnimation(fig=fig, func=update, frames=int(end_time/time_increment), interval=int(time_increment))

ax.axis('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.grid()
plt.show()

plt.plot(np.arange(0, end_time, time_increment), displaced_shape_points[0][:, 1, 1])

plt.show()