classdef RigidBody_Cuboid < RigidBody
    properties
        a, b, c;
    end
    methods
        function obj = RigidBody_Cuboid (rho, a, b, c)
            obj@RigidBody(1,eye(3));
            m = rho*a*b*c;
            Jx = (1/12)*m*(b^2+c^2);
            Jy = (1/12)*m*(c^2+a^2);
            Jz = (1/12)*m*(a^2+b^2);
            J = diag([Jx, Jy, Jz]);
            obj = obj.mass_and_inertia_matrix(m, J);
            obj.density = rho;
            obj.a = a;
            obj.b = b;
            obj.c = c;
        end
        
        function draw(obj, pos, q)
            arguments
                obj;
                pos = zeros(3,1);
                q = [1; 0; 0; 0];
            end
            persistent x000 x001 x011 x010 x100 x101 x111 x110 
            if isempty(x000)
                a2 = obj.a/2; b2 = obj.b/2; c2 = obj.c/2;
                x000 = [ -a2; -b2; -c2 ];
                x001 = [ -a2; -b2;  c2 ];
                x011 = [ -a2;  b2;  c2 ];
                x010 = [ -a2;  b2; -c2 ];
                x100 = [  a2; -b2; -c2 ];
                x101 = [  a2; -b2;  c2 ];
                x111 = [  a2;  b2;  c2 ];
                x110 = [  a2;  b2; -c2 ];
            end
            obj = obj.quaternion(q);
            rotation_matrix = obj.rotation_matrix;
            y000 = rotation_matrix*x000 + pos;
            y001 = rotation_matrix*x001 + pos;
            y011 = rotation_matrix*x011 + pos;
            y010 = rotation_matrix*x010 + pos;
            y100 = rotation_matrix*x100 + pos;
            y101 = rotation_matrix*x101 + pos;
            y111 = rotation_matrix*x111 + pos;
            y110 = rotation_matrix*x110 + pos;
            hold on;
            y = [ y000, y001, y011, y010, y000 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            y = [ y100, y101, y111, y110, y100 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            y = [ y000, y001, y101, y100, y000 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            y = [ y010, y011, y111, y110, y010 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            y = [ y000, y010, y110, y100, y000 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            y = [ y001, y011, y111, y101, y001 ]; fill3(y(1,:), y(2,:), y(3,:), [0.8, 0.8, 0.8]);
            hold off;
        end
    end
end
