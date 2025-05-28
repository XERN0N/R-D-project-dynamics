from SystemModels import *
from Beams import *
from SystemSolvers import *
from Forces import *
import numpy as np

coordinates = ((0, 0, 0), (0, 0, 0.1), (0, 0, 0.2), (0, 0.1, 0.2))
adjacency_matrix = ((0, 1), (1, 2), (2, 3))
time_increment = 0.01
end_time = 2*np.pi
scaling_factor = 10
playback_rate = 1

model = Beam_Lattice()
model.add_beam_edges(coordinates, adjacency_matrix, EN10220.DN_6_06.value, number_of_elements=1)
model.fix_vertices((0,))

model.add_forces({2: lambda t: [1000*np.sin(2*t), 0, 0, 0, 0, 0]})
# stationary_input_shaping_force(model, np.arange(0, end_time, time_increment), stationary_vertex_IDs=(3,), input_vertex_IDs=(1,))
stationary_input_shaping_force_toeplitz_method(model, np.arange(0, end_time, time_increment), stationary_vertex_IDs=(3,), input_vertex_IDs=(1,), solver_timestep=time_increment)

intital_condition_solver = Static(model).solve()
solver = Newmark(model, intital_condition_solver, end_time, time_increment).solve()

animate_solution(solver, playback_rate=playback_rate, scaling_factor=scaling_factor, color='black')
