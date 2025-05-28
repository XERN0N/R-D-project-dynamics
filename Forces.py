from SystemModels import Beam_Lattice
import numpy as np
import numpy.typing as npt
from scipy.interpolate import make_interp_spline as spline
from typing import Literal, Callable

def continous_force_from_array(forces: npt.ArrayLike, timesteps: npt.ArrayLike, cyclic: bool = True) -> Callable[[float], npt.NDArray]:
    """
    Creates a continous force from an array of forces by linearely interpolation.

    Parameters
    ----------
    forces : array_like
        The forces with shape (6, t) where t is the number of timesteps.
    timesteps : array_like
        The timesteps with shape (t,) corresponding to the forces in the array.
    cyclic : bool, optional
        Whether the force should be cyclic after the last timestep. If false, the function returns a zero-force
        after the last timestep.
        
    Returns
    -------
    Callable[[float], npt.NDArray]
        The continous function that takes time as input and returns the interpolated force.
    """
    forces = np.atleast_2d(forces)
    timesteps = np.atleast_1d(timesteps)
    
    if cyclic:
        def continous_force(time: float) -> npt.NDArray:
            time = time % timesteps[-1]
            return np.array([np.interp(time, timesteps, DOF) for DOF in forces])
    else:
        def continous_force(time: float) -> npt.NDArray:
            if time <= timesteps[-1]:
                return np.array([np.interp(time, timesteps, DOF) for DOF in forces])
            else:
                return np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
    return continous_force

def stationary_input_shaping_force_toeplitz_method(model: Beam_Lattice, timesteps: npt.ArrayLike, stationary_vertex_IDs: npt.ArrayLike, input_vertex_IDs: npt.NDArray, solver_timestep: float) -> None:
    """
    Adds a time dependent force to the model that ensures minimal displacement at given vertices in the model.
    It is done using the block toeplitz matrix for the model, the input vector U with actuators

    Parameters
    -----------
    model : Beam_Lattice
        The model to apply the input shaping to.
    timesteps : array_like
        The timesteps of the force vectors to analyze.
    stationary_vertex_IDs : array_like
        The vertex ID(s) of the model to be stationary.
        All the DOF's in the vertex ID(s) will be stationary.
    input_vertex_IDs : array_like
        The vertex ID(s) where the input forces will be applied.
        """
    #error handling ensuring correct input
    length_timesteps = len(timesteps)
    IDs_to_DOFs = lambda IDs: np.repeat(IDs - [np.sum(np.nonzero(model.graph.vs['fixed'])) <i for i in IDs], 6) * 6 + np.tile(np.arange(6), len(IDs))

    #Get vertex IDs for all relevant inputs
    stationary_vertex_IDs = np.atleast_1d(stationary_vertex_IDs)
    input_vertex_IDs = np.atleast_1d(input_vertex_IDs)
    force_vertex_IDs = np.argwhere(np.array(model.graph.vs['force']) != None).reshape(-1)
    input_force_vertex_IDs = np.union1d(input_vertex_IDs, force_vertex_IDs)
    co_location_vertex_IDs = np.intersect1d(input_vertex_IDs, force_vertex_IDs)
    unique_input_vertex_IDs = np.setdiff1d(input_vertex_IDs, force_vertex_IDs)
    unique_force_vertex_IDs = np.setdiff1d(force_vertex_IDs, input_vertex_IDs)

    #Get DOFs from vertex IDs
    stationary_DOFs = IDs_to_DOFs(stationary_vertex_IDs)
    input_DOFs = IDs_to_DOFs(input_vertex_IDs)
    force_DOFs = IDs_to_DOFs(force_vertex_IDs)
    input_forces_DOFs = np.union1d(input_DOFs, force_DOFs)
    co_location_DOFs = np.intersect1d(input_DOFs, force_DOFs)
    unique_input_DOFs = np.setdiff1d(input_DOFs, force_DOFs)
    unique_force_DOFs = np.setdiff1d(force_DOFs, input_DOFs)

    #Get forces from DOFs
    forces = np.empty((length_timesteps, len(force_DOFs))) #2D empty array not tensor

    for i, time in enumerate(timesteps):
        forces[i, :] = model.get_force_vector(time=time)[force_DOFs]
    #forces.flatten() #flatten to fit StateSpace
    
    #Get correct indices to index toeplitz matrix
    unique_force_indices = np.tile(unique_force_DOFs, length_timesteps) + len(input_forces_DOFs)*np.repeat(np.arange(length_timesteps), len(unique_force_DOFs))
    unique_input_indices = np.tile(unique_input_DOFs, length_timesteps) + len(input_forces_DOFs)*np.repeat(np.arange(length_timesteps), len(unique_input_DOFs))
    
    #Get toeplitz matrix and submatrices
    complete_toeplitz = model.get_toeplitz_matrix('receptence', input_DOFs=input_forces_DOFs, output_DOFs=stationary_DOFs, timestep=solver_timestep, number_of_timesteps=length_timesteps)
    unique_input_toeplitz = complete_toeplitz[:, unique_input_indices]
    unique_force_toeplitz = complete_toeplitz[:, unique_force_indices]

    #Solve for the responses unsing r = h1^-1*(-h2*u)
    input_shaped_responses = -np.linalg.pinv(unique_input_toeplitz) @ unique_force_toeplitz @ forces.flatten()

    #Return forces to IDs
    for i, unique_input_vertex_ID in enumerate(unique_input_vertex_IDs):
        unique_input_identifier = np.arange(i*6, i*6+6)
        unique_input_response_indices = np.tile(unique_input_identifier, length_timesteps) + 6*np.repeat(np.arange(length_timesteps), 6)
        unique_input_response = input_shaped_responses[unique_input_response_indices].reshape(length_timesteps, 6).T    
        input_shaped_response_function = continous_force_from_array(unique_input_response, timesteps, cyclic=True)
        model.add_forces({unique_input_vertex_ID: input_shaped_response_function})
        
    if co_location_vertex_IDs is not None:
        for i, co_location_vertex_ID in enumerate(co_location_vertex_IDs): 
            co_location_in_forces = np.where(co_location_vertex_ID == force_vertex_IDs)       
            co_location_force = forces[:, co_location_in_forces*6:co_location_in_forces*+6].T     
            #model.add_forces({co_location_vertex_ID: })
    


def stationary_input_shaping_force(model: Beam_Lattice, timesteps: npt.ArrayLike, stationary_vertex_IDs: npt.ArrayLike, input_vertex_IDs: npt.NDArray, cyclic: bool = True) -> None:
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
    cyclic : bool
        Whether the force should be cyclic after the last timestep. If false, the function returns a zero-force
        after the last timestep.
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
    
    # Adds the response functions to the model.
    for i, input_vertex_ID in enumerate(input_vertex_IDs):
        vertex_response = responses.T[i*6:i*6+6]
        response_function = continous_force_from_array(vertex_response, timesteps, cyclic)
        model.add_forces({input_vertex_ID: response_function})