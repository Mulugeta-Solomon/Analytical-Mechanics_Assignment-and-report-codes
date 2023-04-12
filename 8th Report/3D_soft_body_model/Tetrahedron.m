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

    methods
        function obj = Tetrahedron(i, j, k, l, pi, pj, pk, pl)
            obj.Vertices = [ i, j, k, l ];
            obj.Faces = [ i,j,k; j,l,k; l,i,k; i,l,j ];
            obj.Volume = det( [ pj-pi, pk-pi, pl-pi ] )/6;
            
            ajkl = (pj(2)*pk(3)-pk(2)*pj(3)) + (pk(2)*pl(3)-pl(2)*pk(3)) + (pl(2)*pj(3)-pj(2)*pl(3));
            akli = (pk(2)*pl(3)-pl(2)*pk(3)) + (pl(2)*pi(3)-pi(2)*pl(3)) + (pi(2)*pk(3)-pk(2)*pi(3));
            alij = (pl(2)*pi(3)-pi(2)*pl(3)) + (pi(2)*pj(3)-pj(2)*pi(3)) + (pj(2)*pl(3)-pl(2)*pj(3));
            aijk = (pi(2)*pj(3)-pj(2)*pi(3)) + (pj(2)*pk(3)-pk(2)*pj(3)) + (pk(2)*pi(3)-pi(2)*pk(3));
            
            bjkl = (pj(3)*pk(1)-pk(3)*pj(1)) + (pk(3)*pl(1)-pl(3)*pk(1)) + (pl(3)*pj(1)-pj(3)*pl(1));
            bkli = (pk(3)*pl(1)-pl(3)*pk(1)) + (pl(3)*pi(1)-pi(3)*pl(1)) + (pi(3)*pk(1)-pk(3)*pi(1));
            blij = (pl(3)*pi(1)-pi(3)*pl(1)) + (pi(3)*pj(1)-pj(3)*pi(1)) + (pj(3)*pl(1)-pl(3)*pj(1));
            bijk = (pi(3)*pj(1)-pj(3)*pi(1)) + (pj(3)*pk(1)-pk(3)*pj(1)) + (pk(3)*pi(1)-pi(3)*pk(1));
            
            cjkl = (pj(1)*pk(2)-pk(1)*pj(2)) + (pk(1)*pl(2)-pl(1)*pk(2)) + (pl(1)*pj(2)-pj(1)*pl(2));
            ckli = (pk(1)*pl(2)-pl(1)*pk(2)) + (pl(1)*pi(2)-pi(1)*pl(2)) + (pi(1)*pk(2)-pk(1)*pi(2));
            clij = (pl(1)*pi(2)-pi(1)*pl(2)) + (pi(1)*pj(2)-pj(1)*pi(2)) + (pj(1)*pl(2)-pl(1)*pj(2));
            cijk = (pi(1)*pj(2)-pj(1)*pi(2)) + (pj(1)*pk(2)-pk(1)*pj(2)) + (pk(1)*pi(2)-pi(1)*pk(2));
            
            obj.vector_a = ( 1/(6*obj.Volume))*[ -ajkl; akli; -alij; aijk ];
            obj.vector_b = ( 1/(6*obj.Volume))*[ -bjkl; bkli; -blij; bijk ];
            obj.vector_c = ( 1/(6*obj.Volume))*[ -cjkl; ckli; -clij; cijk ];
            
            a = obj.vector_a;
            b = obj.vector_b;
            c = obj.vector_c;
            
            l = [ a*a', a*b', a*c';
                  b*a', b*b', b*c';
                  c*a', c*b', c*c' ];
            m = [ 2*a*a'+b*b'+c*c', b*a', c*a';
                  a*b', 2*b*b'+c*c'+a*a', c*b';
                  a*c', b*c', 2*c*c'+a*a'+b*b' ];
            
            enum = [ 1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12 ];
            obj.Partial_J_lambda = l(enum, enum);
            obj.Partial_J_mu     = m(enum, enum);
        end

        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l; obj.mu = m;
        end

        function obj = partial_derivaties(obj, ui, uj, uk, ul)
            gamma_u =  [ ui(1); uj(1); uk(1); ul(1) ];
            gamma_v =  [ ui(2); uj(2); uk(2); ul(2) ];
            gamma_w =  [ ui(3); uj(3); uk(3); ul(3) ];
            obj.u_x = obj.vector_a' * gamma_u;
            obj.u_y = obj.vector_b' * gamma_u;
            obj.u_z = obj.vector_c' * gamma_u;
            obj.v_x = obj.vector_a' * gamma_v;
            obj.v_y = obj.vector_b' * gamma_v;
            obj.v_z = obj.vector_c' * gamma_v;
            obj.w_x = obj.vector_a' * gamma_w;
            obj.w_y = obj.vector_b' * gamma_w;
            obj.w_z = obj.vector_c' * gamma_w;
        end

        function obj = calculate_Cauchy_strain(obj, ui, uj, uk, ul)
            obj = obj.partial_derivaties (ui, uj, uk, ul);
            obj.Cauchy_strain = [
                obj.u_x; obj.v_y;  obj.w_z;
                obj.v_z + obj.w_y;
                obj.w_x + obj.u_z;
                obj.u_y + obj.v_x ];
        end

        function energy = partial_strain_potential_energy(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            obj = obj.calculate_Cauchy_strain(ui, uj, uk, ul);
            energy = (1/2) * (obj.Volume) * ...
                ( obj.lambda * ( obj.Cauchy_strain(1) + obj.Cauchy_strain(2) + obj.Cauchy_strain(3) )^2 + ...
                  obj.mu * ( 2*obj.Cauchy_strain(1)^2 + 2*obj.Cauchy_strain(2)^2 + 2*obj.Cauchy_strain(3)^2 + ...
                               obj.Cauchy_strain(4)^2 +   obj.Cauchy_strain(5)^2 +   obj.Cauchy_strain(6)^2 ) );
        end

        function obj = calculate_Green_strain(obj, ui, uj, uk, ul)
            obj = obj.partial_derivaties (ui, uj, uk, ul);
            obj.Green_strain = [ ...
                obj.u_x + (1/2)*(obj.u_x^2 + obj.v_x^2 + obj.w_x^2); ...
                obj.v_y + (1/2)*(obj.u_y^2 + obj.v_y^2 + obj.w_y^2); ...
                obj.w_z + (1/2)*(obj.u_z^2 + obj.v_z^2 + obj.w_z^2); ...
                obj.v_z + obj.w_y + (obj.u_y*obj.u_z + obj.v_y*obj.v_z + obj.w_y*obj.w_z); ...
                obj.w_x + obj.u_z + (obj.u_z*obj.u_x + obj.v_z*obj.v_x + obj.w_z*obj.w_x); ...
                obj.u_y + obj.v_x + (obj.u_x*obj.u_y + obj.v_x*obj.v_y + obj.w_x*obj.w_y) ];
        end

        function energy = partial_strain_potential_energy_Green_strain(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            obj = obj.calculate_Green_strain(ui, uj, uk, ul);
            energy = (1/2) * (obj.Volume) * ...
                ( obj.lambda * ( obj.Green_strain(1) + obj.Green_strain(2) + obj.Green_strain(3) )^2 + ...
                  obj.mu * ( 2*obj.Green_strain(1)^2 + 2*obj.Green_strain(2)^2 + 2*obj.Green_strain(3)^2 + ...
                               obj.Green_strain(4)^2 +   obj.Green_strain(5)^2 +   obj.Green_strain(6)^2 ) );
        end

        function energy = partial_gravitational_potential_energy(obj, disps, grav)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            energy = - (obj.Volume) * grav' *(ui+uj+uk+ul)/4;
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
            mass = obj.Density * obj.Volume;
            I = eye(3);
            obj.Partial_Inertia_Matrix = (mass/20) * ...
                [2*I, I, I, I; I, 2*I, I, I; I, I, 2*I, I; I, I, I, 2*I];
        end

        function [obj, M_p] = partial_inertia_matrix(obj)
            if isempty( obj.Partial_Inertia_Matrix )
                obj = obj.calculate_partial_inertia_matrix;
            end
            M_p = obj.Partial_Inertia_Matrix;
        end

        function obj = calculate_partial_gravitational_vector(obj, g)
            mass = obj.Density * obj.Volume;
            obj.Partial_Gravitational_Vector = mass/4 * [ g; g; g; g ]; 
        end

        function [obj, grav_p] = partial_gravitational_vector(obj, g)
            if isempty( obj.Partial_Gravitational_Vector )
                obj = obj.calculate_partial_gravitational_vector(g);
            end
            grav_p = obj.Partial_Gravitational_Vector;
        end

        function [fi, fj, fk, fl] = nodal_forces_Cauchy_strain(obj, ui, uj, uk, ul)
            a = obj.vector_a;
            b = obj.vector_b;
            c = obj.vector_c;
            obj = obj.calculate_Cauchy_strain (ui, uj, uk, ul);
            mat = obj.lambda*[ones(3,3), zeros(3,3); zeros(3,3), zeros(3,3)] + ...
                      obj.mu*diag([2, 2, 2, 1, 1, 1]);
            Up_e = mat*obj.Cauchy_strain*obj.Volume;
            Up_gammau = [ a, zeros(4,1), zeros(4,1), zeros(4,1), c, b ]*Up_e;
            Up_gammav = [ zeros(4,1), b, zeros(4,1), c, zeros(4,1), a ]*Up_e;
            Up_gammaw = [ zeros(4,1), zeros(4,1), c, b, a, zeros(4,1) ]*Up_e;
            fi = - [ Up_gammau(1); Up_gammav(1); Up_gammaw(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2); Up_gammaw(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3); Up_gammaw(3) ];
            fl = - [ Up_gammau(4); Up_gammav(4); Up_gammaw(4) ];
        end