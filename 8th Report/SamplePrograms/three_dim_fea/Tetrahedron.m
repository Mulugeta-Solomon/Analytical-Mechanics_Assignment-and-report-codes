classdef Tetrahedron
    properties
        Vertices;
        Faces;
        Volume;
        Density; lambda; mu; lambda_vis; mu_vis;
        color;
        vector_a; vector_b; vector_c;
        u_x; u_y; u_z; v_x; v_y; v_z; w_x; w_y; w_z;
        Cauchy_strain;
        Green_strain;
        Partial_J_lambda; Partial_J_mu;
        Partial_Stiffness_Matrix;
        Partial_Damping_Matrix;
        Partial_Inertia_Matrix;
        Partial_Gravitational_Vector;
    end