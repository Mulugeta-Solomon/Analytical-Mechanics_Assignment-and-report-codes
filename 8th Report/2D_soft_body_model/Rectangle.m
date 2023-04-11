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
    methods
        function obj = Rectangle (i, j, k, l, lx, ly, h)
            obj.Vertices = [ i, j, k, l ];
            obj.length_x = lx; obj.length_y = ly;
            obj.Area = lx*ly;
            obj.Thickness = h;
            
            luu = (h/6)*(ly/lx)* ...
                [ 2, -2, -1,  1;
                 -2,  2,  1, -1;
                 -1,  1,  2, -2;
                  1, -1, -2,  2 ];
            luv = (h/4)* ...
                [ 1,  1, -1, -1;
                 -1, -1,  1,  1;
                 -1, -1,  1,  1;
                  1,  1, -1, -1 ];
            lvu = luv';
            lvv = (h/6)*(lx/ly)* ...
                [ 2,  1, -1, -2;
                  1,  2, -2, -1;
                 -1, -2,  2,  1;
                 -2, -1,  1,  2 ];
            muu = 2*luu + lvv; muv = lvu;
            mvu = luv; mvv = 2*lvv + luu;
            
            l = [ luu, luv; lvu, lvv ];
            m = [ muu, muv; mvu, mvv ];
            
            obj.Partial_J_lambda = l([1,5,2,6,3,7,4,8], [1,5,2,6,3,7,4,8]);
            obj.Partial_J_mu     = m([1,5,2,6,3,7,4,8], [1,5,2,6,3,7,4,8]);
        end
        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l; obj.mu = m;
        end

        function obj = calculate_partial_stiffness_matrix(obj)
            obj.Partial_Stiffness_Matrix = ...
                obj.lambda * obj.Partial_J_lambda + obj.mu * obj.Partial_J_mu;
        end

        %function [obj, K_p, Jlambda_p, Jmu_p] = partial_stiffness_matrix(obj)
        function [obj, K_p] = partial_stiffness_matrix(obj)
            if isempty( obj.Partial_Stiffness_Matrix )
                obj = obj.calculate_partial_stiffness_matrix;
            end
            K_p = obj.Partial_Stiffness_Matrix;
            %Jlambda_p = obj.Partial_J_lambda;
            %Jmu_p =     obj.Partial_J_mu;
        end

        function obj = calculate_partial_damping_matrix(obj)
            obj.Partial_Damping_Matrix = ...
                obj.lambda_vis * obj.Partial_J_lambda + obj.mu_vis * obj.Partial_J_mu;
        end

         function [obj, B_p] = partial_damping_matrix(obj)
            if isempty( obj.Partial_Damping_Matrix )
                obj = obj.calculate_partial_damping_matrix;
            end
            B_p = obj.Partial_Damping_Matrix;
        end

        function obj = calculate_partial_inertia_matrix(obj)
            mass = obj.Density * obj.Area * obj.Thickness;
            I = eye(2);
            obj.Partial_Inertia_Matrix = (mass/36) * ...
                [4*I, 2*I, I, 2*I; ...
                 2*I, 4*I, 2*I, I; ...
                 I, 2*I, 4*I, 2*I; ...
                 2*I, I, 2*I, 4*I];
        end

        function [obj, M_p] = partial_inertia_matrix(obj)
            if isempty( obj.Partial_Inertia_Matrix )
                obj = obj.calculate_partial_inertia_matrix;
            end
            M_p = obj.Partial_Inertia_Matrix;
        end

        function energy_density = strain_potential_energy_density(obj, ui, uj, uk, ul, xi, eta)
            lx = obj.length_x; ly = obj.length_y;
            u_x = (1/obj.Area)*( -(ly-eta)*ui(1) +(ly-eta)*uj(1) +eta*uk(1)    -eta *ul(1) );
            u_y = (1/obj.Area)*( -(lx- xi)*ui(1)      -xi *uj(1) +xi *uk(1) +(lx-xi)*ul(1) );
            v_x = (1/obj.Area)*( -(ly-eta)*ui(2) +(ly-eta)*uj(2) +eta*uk(2)    -eta *ul(2) );
            v_y = (1/obj.Area)*( -(lx- xi)*ui(2)      -xi *uj(2) +xi *uk(2) +(lx-xi)*ul(2) );
            energy_density = ...
                obj.lambda * ( u_x + v_y ).^2 + ...
                obj.mu * ( 2*u_x.^2 + 2*v_y.^2 + (u_y + v_x).^2 );
        end