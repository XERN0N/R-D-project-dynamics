from enum import Enum
from typing import TypedDict

class MaterialProperties(TypedDict):
    """
    TypedDict class that holds the material properties of a beam.
    """
    E_modulus: float
    shear_modulus: float
    primary_moment_of_area: float
    secondary_moment_of_area: float
    torsional_constant: float
    density: float
    cross_sectional_area: float

class EN10219(Enum):
    """
    Enum class for material properties of cold formed welded structural 
    hollow sections of non-alloy and fine grain steels (EN10219 beams).
    """
    RHS_30X30X2 = MaterialProperties(
        E_modulus=2.1e11,
        shear_modulus=7.9e10,
        primary_moment_of_area=1.94e-8,
        secondary_moment_of_area=1.02e-8,
        torsional_constant=2.29e-8,
        density=7850,
        cross_sectional_area=1.74e-4
    )

class EN10220(Enum):
    """
    Enum class for material properties of cold formed seamless and welded
    steel tubes (EN 10220 series).
    """
    DN_6_06 = MaterialProperties(
        E_modulus=2.1e11,
        shear_modulus=7.9e10,
        primary_moment_of_area=0.02e-8,
        secondary_moment_of_area=0.02e-8,
        torsional_constant=0.04e-8,
        density=7850,
        cross_sectional_area=0.15e-4
    )