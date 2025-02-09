import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from SystemModels import chain
from SystemSolver import solve_system

def plot_system(
    t: npt.ArrayLike,
    y: npt.ArrayLike,
    plot_position: bool = True,
    plot_velocity: bool = True,
    combined: bool = False
) -> None:
    """
    Plots the position and/or velocity of the system over time.

    Parameters
    ----------
    t : array_like with shape (n,) (1D-vector)
        The time points to plot the system at.
    y : array_like with shape (n, m) (2D-vector)
        The position and velocity of the system.
        Rows are kinematic values, columns are time points.
    plot_position : bool, optional
        Whether to plot the position (default is True).
    plot_velocity : bool, optional
        Whether to plot the velocity (default is True).
    combined : bool, optional
        Whether to combine the position and velocity plots into one plot (default is False).
    
    Returns
    -------
    None
    
    """


    pos, vel = np.split(y,2)
    
    if combined:
        plt.figure()
        if plot_position:
            for i in range(pos.shape[0]):
                plt.plot(t, pos[i], label=f'Position {i+1}')
        if plot_velocity:
            for i in range(vel.shape[0]):
                plt.plot(t, vel[i], label=f'Velocity {i+1}', linestyle='--')
        plt.xlabel('Time (s)')
        plt.ylabel('Position (m) and Velocity (m/s)')
        plt.title('Position and Velocity of the system over time')
        plt.legend()
        plt.grid()
        plt.show()
    else:
        fig, axs = plt.subplots(2, 1, figsize=(10, 8))
        if plot_position:
            for i in range(pos.shape[0]):
                axs[0].plot(t, pos[i], label=f'Position of element {i+1}')
            axs[0].set_xlabel('Time (s)')
            axs[0].set_ylabel('Position (m)')
            axs[0].set_title('Position of the system over time')
            axs[0].legend()
            axs[0].grid()
        if plot_velocity:
            for i in range(vel.shape[0]):
                axs[1].plot(t, vel[i], label=f'Velocity of element {i+1}')
            axs[1].set_xlabel('Time (s)')
            axs[1].set_ylabel('Velocity (m/s)')
            axs[1].set_title('Velocity of the system over time')
            axs[1].legend()
            axs[1].grid()
        plt.tight_layout()
        plt.show()

# Example
if __name__ == '__main__':
    y_0 = np.array([0, 0, 0.1, 0.1, 0, 0,])    
    t = (0, 10)
    Mass, Damping, Stiffness = chain([1, 1, 1,], [3, 2, 1], [0.2, 0.2, 0.2])
    Results_t, Results_y = solve_system(y_0, t, Mass, Damping, Stiffness, 'RK45', 1e-4)
    # Plot the system
    plot_system(Results_t, Results_y, plot_position=True, plot_velocity=True, combined=True)
    
    # Expected output: Two plots showing the position and velocity of the system over time.
    # The position plot should show a sinusoidal curve, and the velocity plot should show a cosine curve. 
    # Both plots should be periodic and have a frequency of 1 Hz.
    # The plots should be displayed in two separate windows.