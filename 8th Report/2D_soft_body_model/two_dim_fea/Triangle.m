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

        function energy = partial_strain_potential_energy(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            obj = obj.calculate_Cauchy_strain(ui, uj, uk);
            energy = (1/2) * (obj.Area * obj.Thickness) * ...
                ( obj.lambda * ( obj.Cauchy_strain(1) + obj.Cauchy_strain(2) )^2 + ...
                  obj.mu * ( 2*obj.Cauchy_strain(1)^2 + 2*obj.Cauchy_strain(2)^2 + obj.Cauchy_strain(3)^2 ) );
        end

         function obj = calculate_Green_strain(obj, ui, uj, uk)
            obj = obj.partial_derivaties (ui, uj, uk);
            obj.Green_strain = [ ...
                obj.u_x + (1/2)*(obj.u_x^2 + obj.v_x^2); ...
                obj.v_y + (1/2)*(obj.u_y^2 + obj.v_y^2); ...
                obj.u_y + obj.v_x + (obj.u_x*obj.u_y + obj.v_x*obj.v_y) ];
        end

        
        function energy = partial_strain_potential_energy_Green_strain(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            obj = obj.calculate_Green_strain(ui, uj, uk);
            energy = (1/2) * (obj.Area * obj.Thickness) * ...
                ( obj.lambda * ( obj.Green_strain(1) + obj.Green_strain(2) )^2 + ...
                  obj.mu * ( 2*obj.Green_strain(1)^2 + 2*obj.Green_strain(2)^2 + obj.Green_strain(3)^2 ) );
        end

         function energy = partial_gravitational_potential_energy(obj, disps, grav)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            energy = - (obj.Density*obj.Area*obj.Thickness) * grav' *(ui+uj+uk)/3;
        end

         function obj = calculate_partial_stiffness_matrix(obj)            
            obj.Partial_Stiffness_Matrix = ...
                obj.lambda * obj.Partial_J_lambda + obj.mu * obj.Partial_J_mu;
        end

        function [obj, K_p] = partial_stiffness_matrix(obj)
            if isempty( obj.Partial_Stiffness_Matrix )
                obj = obj.calculate_partial_stiffness_matrix;
            end
            K_p = obj.Partial_Stiffness_Matrix;
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
            obj.Partial_Inertia_Matrix = (mass/12) * ...
                [2*I, I, I; I, 2*I, I; I, I, 2*I];
        end

         function [obj, M_p] = partial_inertia_matrix(obj)
            if isempty( obj.Partial_Inertia_Matrix )
                obj = obj.calculate_partial_inertia_matrix;
            end
            M_p = obj.Partial_Inertia_Matrix;
        end

          function obj = calculate_partial_gravitational_vector(obj, g)
            mass = obj.Density * obj.Area * obj.Thickness;
            obj.Partial_Gravitational_Vector = mass/3 * [ g; g; g ]; 
        end

         function [fi, fj, fk] = nodal_forces_Cauchy_strain(obj, ui, uj, uk)
            a = obj.vector_a;
            b = obj.vector_b;
            obj = obj.calculate_Cauchy_strain (ui, uj, uk);
            mat = obj.lambda*[1,1,0; 1,1,0; 0,0,0] + obj.mu*[2,0,0; 0,2,0; 0,0,1];
            Up_e = mat*obj.Cauchy_strain*obj.Area*obj.Thickness;
            Up_gammau = [ a, zeros(3,1), b ]*Up_e;
            Up_gammav = [ zeros(3,1), b, a ]*Up_e;
            fi = - [ Up_gammau(1); Up_gammav(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3) ];
        end

           function [fi, fj, fk] = nodal_forces_Green_strain(obj, ui, uj, uk)
            a = obj.vector_a;
            b = obj.vector_b;
            obj = obj.calculate_Green_strain (ui, uj, uk);
            mat = obj.lambda*[1,1,0; 1,1,0; 0,0,0] + obj.mu*[2,0,0; 0,2,0; 0,0,1];
            Up_E = mat*obj.Green_strain*obj.Area*obj.Thickness;
            ux = obj.u_x; uy = obj.u_y; vx = obj.v_x; vy = obj.v_y;
            Up_gammau = [ (1+ux)*a, uy*b, (1+ux)*b+uy*a ]*Up_E;
            Up_gammav = [ vx*a, (1+vy)*b, (1+vy)*a+vx*b ]*Up_E;
            fi = - [ Up_gammau(1); Up_gammav(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3) ];
        end

    end
end