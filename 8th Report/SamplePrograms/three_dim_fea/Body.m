classdef Body
    properties
        numNodalPoints; NodalPoints;
        numTetrahedra; Tetrahedra;
        BoundaryFaces;
        numSubRegions; SubRegions;
        numSurfaces; Surfaces;
        Density; lambda; mu; lambda_vis; mu_vis;
        strain_potential_energy;
        gravitational_potential_energy;
        J_lambda; J_mu;
        Stiffness_Matrix;
        Damping_Matrix;
        Inertia_Matrix;
        Gravitational_Vector;
    end