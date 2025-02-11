import numpy as np
import numpy.typing as npt
from typing import Callable
from scipy.integrate import solve_ivp
from SystemModels import chain

def harmonic_load(t: float) -> float:
    """
    Defines a harmonic load as a function of time.

    Parameters
    ----------
    t : float
        The time point to evaluate the load at.

    Returns
    -------
    float
        The load at the given time point.
    """
    return np.sin(t)

def equation_of_motions(
    t: float,
    y: npt.NDArray, 
    M: npt.NDArray, 
    C: npt.NDArray, 
    K: npt.NDArray,
    B_2: npt.NDArray,
    u: Callable[[float], npt.NDArray],
    ) -> tuple[npt.ArrayLike]:
    """
    Calculates the velocity and acceleration of the system at a given time.
    The system is defined by the mass matrix, stiffness matrix, and damping matrix.

    Parameters
    ----------
    t : array_like with shape (n,) (1D-vector)
        The time points to evaluate the system at.
    y : array_like with shape (n,) (1D-vector)
        The position and velocity of the system.
    M : array_like with shape (n, n) (2D-vector)
        The mass matrix of the system.
    C : array_like with shape (n, n) (2D-vector)
        The damping matrix of the system.
    K : array_like with shape (n, n) (2D-vector)
        The stiffness matrix of the system.
    B_2 : array_like with shape (n, n) (2D-vector)
        The matrix for the external force mapping.
    u : array_like with shape (n,) (1D-vector)
        The external force applied to the system.

    Returns
    -------
    tuple of np arrays with shape (n,) (1D-vector)
        The velocity and acceleration of the system.
    """

    #Splitting the 1D input with [position, velocity] into two 1D arrays
    y, ydot = np.split(y, 2)

    # Equation of motion for the system
    # [M]*y'' + [C]*y' + [K]*y -[B_2]*u = 0
    # Isolate ddy/dt:
    # y'' = [M]^-1 * (-[C]*y' - [K]*y + [B_2]*u)

    ydotdot = np.linalg.inv(M) @ (-C @ ydot - K @ y + B_2 @ u(t))

    # Returning velocity and acceleration as a concatenated 1D array
    return np.concatenate((ydot, ydotdot))

#Example
if __name__ == "__main__":
    y_0 = np.array([0, 0, 1, 0, 0, 0])
    t = (0, 10)
    t_points = np.linspace(*t, 1000)
    Mass, Damping, Stiffness = chain([1, 1, 1], [3, 2, 1], [0.2, 0.2, 0.2])
    load = lambda t: np.array([1]) * np.sin(t*1)
    load_map = np.array([[0, 0, 1]]).T
    Results = solve_ivp(
        fun=equation_of_motions, 
        t_span=t, 
        y0=y_0, 
        args=(Mass, Damping, Stiffness, load_map, load),
        t_eval=t_points
    )
    Results_y = Results.y
    Results_t = Results.t
    combined_results = np.hstack((Results_y[:, :3], Results_y[:, -3:]))
    np.set_printoptions(precision=3)
    print(f"Combined results (6x??): \n{combined_results}")
    print(f"The timepoints were: \n{Results_t}")
    print(Results_y.shape)