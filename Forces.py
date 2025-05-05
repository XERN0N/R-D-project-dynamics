from SystemModels import Beam_Lattice
import numpy as np
import numpy.typing as npt
from typing import Literal

def stationary_input_shaping_force(model: Beam_Lattice, timesteps: npt.ArrayLike, stationary_vertex_IDs: npt.ArrayLike, input_vertex_IDs: npt.NDArray, response_shape: Literal['cyclic', 'non-cyclic'] = 'cyclic') -> None:
    """
    Adds a time dependent force to the model that ensures minimal displacement at given vertices in the model.
    If the system is overdetermined (more outputs than inputs) the input vector which minimizes the L2 norm of
    the output vector is used. If underdetermined (more inputs than outputs) the L2 norm of the response vector
    is minimized.

    Parameters
    ----------
    model : Beam_Lattice
        The model to apply the input shaping to.
    timesteps : array_like
        The timesteps of the force vectors to analyze.
    stationary_vertex_IDs : array_like
        The vertex ID(s) of the model to be stationary.
        All the DOF's in the vertex ID(s) will be stationary.
    input_vertex_IDs : array_like
        The vertex ID(s) where the input forces will be applied.
    response_shape : str
        Extrapolated response type. Can be either 'cyclic' or 'non-cyclic' where non-cyclic adds zero-padding to the response and 'cyclic'
        repeats the response with a period corresponding to the length of the timesteps.
    """    
    stationary_vertex_IDs = np.atleast_1d(stationary_vertex_IDs)
    input_vertex_IDs = np.atleast_1d(input_vertex_IDs)
    IDs_to_DOFs = lambda IDs: np.repeat(IDs - [np.sum(np.nonzero(model.graph.vs['fixed'])) < i for i in IDs], 6) * 6 + np.tile(np.arange(6), len(IDs))

    # Get the force vector from the model for all time steps.
    force_vertex_IDs = np.argwhere(np.array(model.graph.vs['force']) != None).reshape(-1)
    force_DOFs = IDs_to_DOFs(force_vertex_IDs)
    forces = np.empty((len(timesteps), len(force_DOFs), 1))
    for i, time in enumerate(timesteps):
        forces[i, :] = model.get_force_vector(time=time)[force_DOFs].reshape(-1, 1)

    # Perform transformation of force vector.
    forces_in_frequency_domain = np.fft.fft(forces.reshape(len(timesteps), len(force_DOFs)), axis=0).reshape(len(timesteps), len(force_DOFs), 1)
    frequencies = np.fft.fftfreq(len(timesteps), timesteps[-1]/len(timesteps))

    """Solves for the response function in the laplace domain."""
    # Determines which DOF are accociated with the input and stationary vertices.
    input_DOFs = IDs_to_DOFs(input_vertex_IDs)
    stationary_DOFs = IDs_to_DOFs(stationary_vertex_IDs)
    # Gets the transfer matrix and picks the associated transfer functions.
    transfer_matrix_rows = np.tile(np.atleast_2d(stationary_DOFs).T, (len(frequencies), 1, 1))
    combined_forces_DOFs = np.sort(np.concatenate((input_DOFs, force_DOFs)))
    transfer_matrix_columns = np.tile(np.atleast_2d(combined_forces_DOFs), (len(frequencies), 1, 1))
    transfer_function = np.take_along_axis(
                        np.take_along_axis(model.get_transfer_matrix('receptence', frequencies), 
                                                transfer_matrix_rows,    axis=1), 
                                                transfer_matrix_columns, axis=2)
    # Rearanging the transfer matrix such that the last columns relates to the response forces and the first to the forces.
    boolean_mask_of_columns_related_to_response_forces = np.isin(combined_forces_DOFs, input_DOFs)
    transfer_function_column_indices = np.arange(len(input_DOFs) + len(force_DOFs))
    rearanged_column_indices = np.concatenate((transfer_function_column_indices[~boolean_mask_of_columns_related_to_response_forces], 
                                               transfer_function_column_indices[ boolean_mask_of_columns_related_to_response_forces]))
    transfer_function = transfer_function[:, :, rearanged_column_indices]

    # Isolates the response function in the laplace domain.
    response_in_frequency_domain = np.linalg.pinv(transfer_function[:, :, len(force_DOFs):]) @ (-transfer_function[:, :, :len(force_DOFs)] @ forces_in_frequency_domain)

    # Apply inverse z-transform to the response function.
    responses = np.fft.ifft(response_in_frequency_domain.reshape(len(timesteps), len(input_DOFs)), axis=0).real
    
    for i, input_vertex_ID in enumerate(input_vertex_IDs):
        vertex_response = responses.T[i*6:i*6+6]
        match response_shape:
            case 'cyclic':
                def response_function(time: float) -> npt.NDArray:
                    time = time % timesteps[-1]
                    return np.array([np.interp(time, timesteps, DOF_response) for DOF_response in vertex_response])
            case 'non-cyclic':
                def response_function(time: float) -> npt.NDArray:
                    if time <= timesteps[-1]:
                        return np.array([np.interp(time, timesteps, DOF_response) for DOF_response in vertex_response])
                    else:
                        return np.array([0, 0, 0, 0, 0, 0])
                    
        model.add_forces({input_vertex_ID: response_function})