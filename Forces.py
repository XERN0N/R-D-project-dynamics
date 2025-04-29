from SystemModels import Beam_Lattice
import numpy as np
import numpy.typing as npt

def stationary_input_shaping_force(model: Beam_Lattice, timesteps: npt.ArrayLike, stationary_vertex_IDs: npt.ArrayLike, input_vertex_IDs: npt.NDArray) -> None:
    """
    Adds a time dependent force to the model that ensures minimal displacement at given vertices in the model.
    If the system is overdetermined (more outputs than inputs) the input vector which minimizes the L2 norm of
    the output vector is used. Underdetermined systems are not supported.

    Parameters
    ----------
    model : Beam_Lattice
        The model to apply the input shaping to.
    timesteps : array_like
        The timesteps of the force vectors to analyze.
    stationary_vertex_IDs : array_like
        The vertex ID(s) of the model to be stationary where len(stationary_vertex_IDs) >= len(input_vertex_IDs).
        All the DOF's in the vertex ID(s) will be stationary.
    input_vertex_IDs : array_like
        The vertex ID(s) where the input forces will be applied.
    """
    if len(stationary_vertex_IDs) < len(input_vertex_IDs):
        raise ValueError(f"Cannot solve for a underdetermined system (stationary vertices ({len(stationary_vertex_IDs)}) < input vertices ({len(input_vertex_IDs)})).")
    
    # Get the force vector from the model for all time steps.
    pass

    # Get the transfer matrix from the model.
    pass

    # Perform transformation of force vector.
    pass

    # Solve for the input vector than minimizes the L2 norm of the output vector.
    pass

    # Add the new force to the model.
    pass