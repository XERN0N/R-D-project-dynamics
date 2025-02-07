import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from SystemModels import chain

def equation_of_motions(
    t: npt.ArrayLike,
    y: npt.ArrayLike, 
    M: npt.ArrayLike, 
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

def solve_system(
    y_ic: npt.ArrayLike,
    t_duration: npt.ArrayLike,
    M: npt.ArrayLike,
    C: npt.ArrayLike,
    K: npt.ArrayLike,
    sol_method: str = 'RK45',
    rel_tol: float = None,
) -> tuple[npt.ArrayLike, npt.ArrayLike]:
    """
    Solves the system of differential equations defined by the mass, damping, and stiffness matrices using chain function.

    Parameters
    ----------
    y_ic : array_like with shape (n,) (1D-vector)
        The initial conditions of the system.
    t_duration : array_like with shape (2,) (1D-vector)
        The time span for the solution.
    M : array_like with shape (n, n) (2D-vector)
        The mass matrix of the system.
    C : array_like with shape (n, n) (2D-vector)
        The damping matrix of the system.
    K : array_like with shape (n, n) (2D-vector)
        The stiffness matrix of the system.
    sol_method : str, optional
        The integration method to use (default is 'RK45').
    rel_tol : float, optional
        The relative tolerance for the solver (default is None).

    Returns
    -------
    tuple of np arrays
        The time points and the solution of the system.
    """
    # This is a series of if-else statements to handle unspecified sol_method and rel_tol.
    if sol_method and rel_tol is None:
        Ivp_Results = solve_ivp(equation_of_motions, t_duration, y_ic, method=sol_method, args=(M, C, K))
    elif rel_tol is None:
        Ivp_Results = solve_ivp(equation_of_motions, t_duration, y_ic, method=sol_method, args=(M, C, K))
    elif sol_method is None:
        Ivp_Results = solve_ivp(equation_of_motions, t_duration, y_ic, method='RK45', rtol=rel_tol, args=(M, C, K))
    else:
        Ivp_Results = solve_ivp(equation_of_motions, t_duration, y_ic, method=sol_method, rtol=rel_tol, args=(M, C, K))
    
    return (Ivp_Results.t, Ivp_Results.y)



#Example
if __name__ == "__main__":
    y_0 = np.array([0, 0, 1, 0, 0, 0])
    t = (0, 10)
    Mass, Damping, Stiffness = chain([1, 1, 1], [3, 2, 1], [0.2, 0.2, 0.2])
    Results_t, Results_y = solve_system(y_0, t, Mass, Damping, Stiffness, 'RK23', 0.01)
    
    combined_results = np.hstack((Results_y[:, :3], Results_y[:, -3:]))
    np.set_printoptions(precision=3)
    print(f"Combined results (6x??): \n{combined_results}")
    print(f"The timepoints were: \n{Results_t}")
    print(Results_y.shape)