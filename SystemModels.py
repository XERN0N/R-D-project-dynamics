import numpy as np
import numpy.typing as npt

def chain(num_of_masses: int, end_connection: bool = False) -> tuple[npt.NDArray, npt.NDArray]:
    stifness_matrix = np.eye(num_of_masses)



def mass_spring_damper(m: float, k: float, c: float) -> tuple[npt.NDArray, npt.NDArray]:
    A = np.array([[0, 1], [-k/m, -c/m]])
    B = np.array([[0], [1/m]])
    return A, B