import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from SystemModels import chain
from SystemSolver import *
from scipy.fft import fft, fftfreq

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
                plt.plot(t, pos[i], label=f'Position {i+1}', marker='x')
        if plot_velocity:
            for i in range(vel.shape[0]):
                plt.plot(t, vel[i], label=f'Velocity {i+1}', linestyle='--', marker='o')
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
                axs[0].plot(t, pos[i], label=f'Position of element {i+1}', marker='x')
            axs[0].set_xlabel('Time (s)')
            axs[0].set_ylabel('Position (m)')
            axs[0].set_title('Position of the system over time')
            axs[0].legend()
            axs[0].grid()
        if plot_velocity:
            for i in range(vel.shape[0]):
                axs[1].plot(t, vel[i], label=f'Velocity of element {i+1}', marker='o')
            axs[1].set_xlabel('Time (s)')
            axs[1].set_ylabel('Velocity (m/s)')
            axs[1].set_title('Velocity of the system over time')
            axs[1].legend()
            axs[1].grid()
        plt.tight_layout()
        plt.show()

def plot_fft(
    t: npt.ArrayLike,
    y: npt.ArrayLike
) -> None:
    """
    Plots the FFT of the system's position and velocity over time.

    Parameters
    ----------
    t : array_like with shape (n,) (1D-vector)
        The time points of the system.
    y : array_like with shape (n, m) (2D-vector)
        The position and velocity of the system.
    
    Returns
    -------
    None
    """
    N = len(t)
    T = t[1] - t[0]  # Assuming uniform sampling
    yf = fft(y, axis=1)
    xf = fftfreq(N, T)[:N//2]

    plt.figure()
    for i in range(y.shape[0]):
        plt.plot(xf, 2.0/N * np.abs(yf[i, :N//2]), label=f'Component {i+1}')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('FFT of the system')
    plt.legend()
    plt.grid()
    plt.show()

# Example
if __name__ == '__main__':
    y_0 = np.array([0, 0, 0, 0, 0, 0])
    t = (0, 50)
    t_points = np.linspace(*t, 2048)
    Mass, Damping, Stiffness = chain([1, 1, 1], [1, 1, 1], [100, 100, 100])
    load = lambda t: np.array([np.sin(t*1)*1, np.sin(t*1)*-1]) 
    load_map = np.array([[0, 0, 1], [0, 1, 0]]).T
    Results = solve_ivp(
        fun=equation_of_motions, 
        t_span=t, 
        y0=y_0, 
        args=(Mass, Damping, Stiffness, load_map, load),
        t_eval=t_points
    )
    Results_y = Results.y[:3, 20:]
    Results_t = Results.t[20:]

    plot_fft(Results_t, Results_y)  # Plot FFT for all components

    plot_system(Results_t, Results_y, plot_position=True, plot_velocity=True, combined=False)