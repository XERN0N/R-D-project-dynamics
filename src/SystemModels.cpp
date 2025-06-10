#include <cblas.h>
#include <SystemModels.h>
#include <math.h>
#include <algorithm>
using namespace std;

Edge::Edge(Vertex& start_vertex, Vertex& end_vertex, double E_modulus, 
           double shear_modulus, double primary_moment_of_area, double secondary_moment_of_area, 
           double torsional_constant, double density, double cross_sectional_area)
           : start_vertex(start_vertex), end_vertex(end_vertex) {
    
    // Calculates the edge length.
    double edge_vector[3] = {0, 0, 0};
    cblas_daxpy(3, 1, start_vertex.coordinate, 1, edge_vector, 1);
    cblas_daxpy(3, -1, end_vertex.coordinate, 1, edge_vector, 1);
    edge_length = cblas_dnrm2(3, edge_vector, 1);

    // Short handing the beam properties.
    double& E = E_modulus; double& G = shear_modulus; double& I_z = primary_moment_of_area; double& L = edge_length;
    double& I_y = secondary_moment_of_area; double& J = torsional_constant; double& rho = density; double& A = cross_sectional_area;
    double I0 = I_y + I_z;
    
    // Defines the mass matrix.
    double temporary_edge_mass_matrix[] = {140,     0,     0,        0,            0,            0,  70,     0,     0,        0,            0,            0,
                                             0,   156,     0,        0,            0,         22*L,   0,    54,     0,        0,            0,        -13*L,
                                             0,     0,   156,        0,        -22*L,            0,   0,     0,    54,        0,         13*L,            0,
                                             0,     0,     0, 140*I0/A,            0,            0,   0,     0,     0,  70*I0/A,            0,            0,
                                             0,     0, -22*L,        0,  4*pow(L, 2),            0,   0,     0, -13*L,        0, -3*pow(L, 2),            0,
                                             0,  22*L,     0,        0,            0,  4*pow(L, 2),   0,  13*L,     0,        0,            0, -3*pow(L, 2),
                                            70,     0,     0,        0,            0,            0, 140,     0,     0,        0,            0,            0,
                                             0,    54,     0,        0,            0,         13*L,   0,   156,     0,        0,            0,        -22*L,
                                             0,     0,    54,        0,        -13*L,            0,   0,     0,   156,        0,         22*L,            0,
                                             0,     0,     0,  70*I0/A,            0,            0,   0,     0,     0, 140*I0/A,            0,            0,
                                             0,     0,  13*L,        0, -3*pow(L, 2),            0,   0,     0,  22*L,        0,  4*pow(L, 2),            0,
                                             0, -13*L,     0,        0,            0, -3*pow(L, 2),   0, -22*L,     0,        0,            0,  4*pow(L, 2)};
    cblas_dscal(144, rho*A*L/420, temporary_edge_mass_matrix, 1); 
    copy(temporary_edge_mass_matrix, temporary_edge_mass_matrix + 144, edge_mass_matrix);

    // Defines the stiffness matrix.
    double temporary_element_stiffness_matrix[] = { E*A/L,                   0,                   0,      0,                  0,                  0, -E*A/L,                   0,                   0,      0,                  0,                  0,
                                                        0,  12*E*I_z/pow(L, 3),                   0,      0,                  0,  6*E*I_z/pow(L, 2),      0, -12*E*I_z/pow(L, 3),                   0,      0,                  0,  6*E*I_z/pow(L, 2),
                                                        0,                   0,  12*E*I_y/pow(L, 3),      0, -6*E*I_y/pow(L, 2),                  0,      0,                   0, -12*E*I_y/pow(L, 3),      0, -6*E*I_y/pow(L, 2),                  0,
                                                        0,                   0,                   0,  G*J/L,                  0,                  0,      0,                   0,                   0, -G*J/L,                  0,                  0,
                                                        0,                   0,  -6*E*I_y/pow(L, 2),      0,          4*E*I_y/L,                  0,      0,                   0,   pow(6*E*I_y/L, 2),      0,          2*E*I_y/L,                  0,
                                                        0,   6*E*I_z/pow(L, 2),                   0,      0,                  0,          4*E*I_z/L,      0,  -6*E*I_z/pow(L, 2),                   0,      0,                  0,          2*E*I_z/L,
                                                   -E*A/L,                   0,                   0,      0,                  0,                  0,  E*A/L,                   0,                   0,      0,                  0,                  0,
                                                        0, -12*E*I_z/pow(L, 3),                   0,      0,                  0, -6*E*I_z/pow(L, 2),      0,  12*E*I_z/pow(L, 3),                   0,      0,                  0, -6*E*I_z/pow(L, 2),
                                                        0,                   0, -12*E*I_y/pow(L, 3),      0,  6*E*I_y/pow(L, 2),                  0,      0,                   0,  12*E*I_y/pow(L, 3),      0,  6*E*I_y/pow(L, 2),                  0,
                                                        0,                   0,                   0, -G*J/L,                  0,                  0,      0,                   0,                   0,  G*J/L,                  0,                  0,
                                                        0,                   0,  -6*E*I_y/pow(L, 2),      0,          2*E*I_y/L,                  0,      0,                   0,   6*E*I_y/pow(L, 2),      0,          4*E*I_y/L,                  0,
                                                        0,   6*E*I_z/pow(L, 2),                   0,      0,                  0,          2*E*I_z/L,      0,  -6*E*I_z/pow(L, 2),                   0,      0,                  0,          4*E*I_z/L};
    
    copy(temporary_element_stiffness_matrix, temporary_element_stiffness_matrix + 144, edge_stiffness_matrix);
}
