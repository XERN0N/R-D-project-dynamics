#ifndef SYSTEMMODELS_H
#define SYSTEMMODELS_H

struct Vertex {
    double coordinate[3];
    bool fixed = false;
    double (*force_function)(double) = nullptr;
};

class Edge {
public:
    double edge_mass_matrix[144];
    double edge_stiffness_matrix[144];
    Vertex& start_vertex;
    Vertex& end_vertex;
    double edge_length;

    Edge(Vertex& start_vertex, Vertex& end_vertex, double E_modulus, 
         double shear_modulus, double primary_moment_of_area, double secondary_moment_of_area, 
         double torsional_constant, double density, double cross_sectional_area);
};

#endif