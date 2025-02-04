import numpy as np
import numpy.typing as npt

def chain(num_of_masses: int, end_connection: bool = False) -> tuple[npt.NDArray, npt.NDArray]:
    stifness_matrix = np.eye(num_of_masses)