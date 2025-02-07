import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from SystemModels import chain

def equation_of_motions(
    t: float,
    y: npt.ArrayLike, 
    M: npt,ArrayLike, 
    C: npt.ArrayLike, 
    K: npt.ArrayLike
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

    Returns
    -------
    tuple of np arrays with shape (n,) (1D-vector)
        The velocity and acceleration of the system.
    """
    #Splitting the 1D input with [position, velocity] into two 1D arrays
    y, ydot = np.split(y, 2)

    # Equation of motion for the system
    # [M]*y'' + [C]*y' + [K]*y = 0
    # Isolate ddy/dt:
    # y'' = [M]^-1 * (-[C]*y' - [K]*y)

    ydotdot = np.linalg.inv(M) @ (-C @ydot - K @ y)

    # Returning velocity and acceleration as a concatenated 1D array
    return np.concatenate((ydot, ydotdot))


y_ic = np.array([0, 0, 0, 0, 0, 1])
t_duration = (10)
M, C, K = chain([1, 1, 1], [3, 2, 1], [0.2, 0.2, 0.2])
Ivp_Results = solve_ivp(equation_of_motions, t_duration, y_ic, method='RK45')

print(Ivp_Results.y)
#Example
