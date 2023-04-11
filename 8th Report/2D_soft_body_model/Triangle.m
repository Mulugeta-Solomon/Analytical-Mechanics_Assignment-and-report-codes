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
    methods
        function obj = Triangle(i, j, k, pi, pj, pk, h)
            obj.Vertices = [ i, j, k ];
            obj.Area = det( [ pj-pi, pk-pi ] )/2;
            obj.Thickness = h;
            obj.vector_a = ( 1/(2*obj.Area))*[ pj(2)-pk(2); pk(2)-pi(2); pi(2)-pj(2) ];
            obj.vector_b = (-1/(2*obj.Area))*[ pj(1)-pk(1); pk(1)-pi(1); pi(1)-pj(1) ];
            
            a = obj.vector_a;
            b = obj.vector_b;
            vol = obj.Area * obj.Thickness;
            luu = vol*a*a'; luv = vol*a*b';
            lvu = vol*b*a'; lvv = vol*b*b';
            muu = 2*luu + lvv; muv = lvu;
            mvu = luv; mvv = 2*lvv + luu;
            
            l = [ luu, luv; lvu, lvv ];
            m = [ muu, muv; mvu, mvv ];
    
            obj.Partial_J_lambda = l([1,4,2,5,3,6], [1,4,2,5,3,6]);
            obj.Partial_J_mu     = m([1,4,2,5,3,6], [1,4,2,5,3,6]);
        end
        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l; obj.mu = m;
        end

        function obj = partial_derivaties(obj, ui, uj, uk)
            gamma_u =  [ ui(1); uj(1); uk(1) ];
            gamma_v =  [ ui(2); uj(2); uk(2) ];
            obj.u_x = obj.vector_a' * gamma_u;
            obj.u_y = obj.vector_b' * gamma_u;
            obj.v_x = obj.vector_a' * gamma_v;
            obj.v_y = obj.vector_b' * gamma_v;
        end

        function obj = calculate_Cauchy_strain(obj, ui, uj, uk)
            obj = obj.partial_derivaties (ui, uj, uk);
            obj.Cauchy_strain = [ obj.u_x; obj.v_y;  obj.u_y + obj.v_x ];
        end
