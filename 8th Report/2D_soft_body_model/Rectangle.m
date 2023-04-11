classdef Rectangle
    properties
        Vertices;
        length_x; length_y;
        Area;
        Thickness;
        Density; lambda; mu; lambda_vis; mu_vis;
        color;
        Partial_J_lambda; Partial_J_mu;
        Partial_Stiffness_Matrix;
        Partial_Damping_Matrix;
        Partial_Inertia_Matrix;
    end