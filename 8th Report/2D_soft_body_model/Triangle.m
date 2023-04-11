classdef Triangle
    properties
        Vertices;
        Area;
        Thickness;
        Density; lambda; mu; lambda_vis; mu_vis;
        color;
        vector_a; vector_b;
        u_x; u_y; v_x; v_y;
        Cauchy_strain;
        Green_strain;
        Partial_J_lambda; Partial_J_mu;
        Partial_Stiffness_Matrix;
        Partial_Damping_Matrix;
        Partial_Inertia_Matrix;
        Partial_Gravitational_Vector;
    end
    
