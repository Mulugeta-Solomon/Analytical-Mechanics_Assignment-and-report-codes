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