classdef RigidBody
    properties
        density;
        mass;
        inertia_matrix;
        inertia_matrix_inverse;
        rotation_matrix;
        omega;
        q;
        dotq;
        H;
    end
    methods
        function obj = RigidBody (m, J)
            obj.mass = m;
            obj.inertia_matrix = J;
            obj.inertia_matrix_inverse = inv(J);
        end
        
        function obj = mass_and_inertia_matrix(obj, m, J)
            obj.mass = m;
            obj.inertia_matrix = J;
            obj.inertia_matrix_inverse = inv(J);
        end
        
        function obj = quaternion (obj, q)
            obj.q = q;
            q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
            obj.rotation_matrix = ...
                [ 2*(q0^2+q1^2)-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2); ...
                  2*(q1*q2+q0*q3), 2*(q0^2+q2^2)-1, 2*(q2*q3-q0*q1); ...
                  2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(q0^2+q3^2)-1 ];
            obj.H = ...
                [ -q1,  q0,  q3, -q2; ...
                  -q2, -q3,  q0,  q1; ...
                  -q3,  q2, -q1,  q0 ];
        end
        
        function obj = dot_quaternion(obj, dotq)
            obj.dotq = dotq;
            obj.omega = 2*obj.H*dotq;
        end
        
        function ddotq = calculate_ddotq (obj, q, dotq, alpha, external_torque)
            obj = obj.quaternion(q);
            obj = obj.dot_quaternion(dotq);
            r = dotq'*dotq + 2*alpha*q'*dotq + (1/2)*alpha^2*(q'*q-1);
            Hdotq = obj.H*dotq;
            JHdotq = obj.inertia_matrix*Hdotq;
            ddotq = -r*q - 2*(obj.H)'*obj.inertia_matrix_inverse*(cross(Hdotq,JHdotq)-(1/4)*external_torque);
        end
    end
end
