import numpy as np
from scipy.integrate import solve_ivp
from SystemModels import chain

def equation_of_motions(t, y, M, K, C) -> None:
    position, velocity = np.split(y, 2)

    acceleration = np.linalg.inv(M) @ (C @ velocity - K @ position)

    return np.concatenate((velocity, acceleration))

initial_condition = np.array([0, 0, 0, 0, 0, 1])
t_span = (0, 1)
M, K, C = chain([1, 1, 1], [3, 2, 1], [5, 4, 3])

results = solve_ivp(equation_of_motions, t_span, initial_condition, args=(M, K, C))

velocity, acceleration = np.split(results.y[:, -1], 2)

print(velocity)