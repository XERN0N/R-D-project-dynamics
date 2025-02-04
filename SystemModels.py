import numpy as np
import numpy.typing as npt

def chain(masses: npt.ArrayLike, springs: npt.ArrayLike, dampers: npt.ArrayLike = None) -> tuple[npt.NDArray, npt.NDArray] | tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Creates the mass and spring matrix with size n and also the damper matrix if supplied.

    Parameters
    ----------
    masses : array_like with shape (n,)
        The mass of the lumped masses in the system.
    springs : array_like with shape (n,) | (n+1,)
        The spring constants between the lumped masses.
    dampers : array_like with shape (n,) | (n+1,), optional
        The damper constants of the dampers in between. Ignored by default.

    Returns
    -------
    tuple of np arrays with shape (n, n)
        Contains the mass matrix, spring matrix, and damper matrix if provided.
    """
    springs = np.asarray(springs)
    matrix_size = len(masses)
    
    # Creating mass matrix.
    mass_matrix = np.eye(matrix_size)
    np.fill_diagonal(mass_matrix, masses)
    
    def create_matrix(actuator: npt.NDArray) -> npt.NDArray:
        matrix = np.zeros((matrix_size, matrix_size))
        diag_indicies = np.diag_indices_from(matrix)
        # If there is no actuator on the end.
        if matrix_size == len(actuator):
            # Main diagonal.
            matrix[diag_indicies] = actuator + np.concatenate((actuator[1:], [0]))
            # Sub diagonals.
            matrix[1:, :-1][np.diag_indices(matrix_size - 1)] = -actuator[1:]
            matrix[:-1, 1:][np.diag_indices(matrix_size - 1)] = -actuator[1:]
        # If there is an actuator on the end.
        elif matrix_size + 1 == len(actuator):
            # Main diagonals.
            matrix[diag_indicies] = actuator[:-1] + actuator[1:]
            # Sub diagonals.
            matrix[1:, :-1][np.diag_indices(matrix_size - 1)] = -actuator[1:-1]
            matrix[:-1, 1:][np.diag_indices(matrix_size - 1)] = -actuator[1:-1]
        else:
            raise ValueError(f"Length of mass and spring vectors aren't compatible (len(masses) = {matrix_size}, len(springs) = {len(actuator)}).")
        return matrix
        
    # Creating the stiffness matrix.
    stiffness_matrix = create_matrix(springs)
    if dampers is not None:
        # Creating the damper matrix.
        dampers = np.asarray(dampers)
        damper_matrix = create_matrix(dampers)
        return mass_matrix, stiffness_matrix, damper_matrix
    
    return mass_matrix, stiffness_matrix

# Example.
if __name__ == "__main__":
    masses = [1, 2, 3]
    springs = [3, 5, 2]
    dampers = [5, 3, 7, 8]
    matrices = chain(masses, springs, dampers)
    print(f"Mass matrix:\n{matrices[0]}\n")
    print(f"Spring matrix:\n{matrices[1]}\n")
    print(f"Damper matrix:\n{matrices[2]}")