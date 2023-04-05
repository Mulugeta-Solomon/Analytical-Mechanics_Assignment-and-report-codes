classdef Open_Mechanism_Two_DOF
    properties
        link1;
        link2;
        base_position;
        gravity;
        theta1; theta2;
        omega1; omega2;
        C1; S1; C2; S2; C12; S12;
    end
    methods
        function obj = Open_Mechanism_Two_DOF (l1, l2, pos, grav)
            obj.link1 = l1;
            obj.link2 = l2;
            obj.base_position = pos;
            obj.gravity = grav;
        end
        
        function obj = joint_angles (obj, theta, omega)
            obj.theta1 = theta(1); obj.theta2 = theta(2);
            obj.C1 = cos(obj.theta1); obj.S1 = sin(obj.theta1);
            obj.C2 = cos(obj.theta2); obj.S2 = sin(obj.theta2);
            obj.C12 = obj.C1*obj.C2 - obj.S1*obj.S2; % cos(theta1+theta2);
            obj.S12 = obj.S1*obj.C2 + obj.C1*obj.S2; % sin(theta1+theta2);
            obj.omega1 = omega(1); obj.omega2 = omega(2);
        end
        
        function draw (obj)
            O = obj.base_position;
            P1 = obj.base_position + obj.link1.length*[ obj.C1; obj.S1 ];
            P2 = P1 + obj.link2.length*[ obj.C12; obj.S12 ];
            x = [O(1); P1(1); P2(1)];
            y = [O(2); P1(2); P2(2)];
            plot(x,y,'k-');
        end

        function tip = tip_point (obj)
            tip = obj.base_position + ...
                  obj.link1.length*[obj.C1; obj.S1] + ...
                  obj.link2.length*[obj.C12; obj.S12];
        end
        
        function J = Jacobian (obj)
            l1 = obj.link1.length;
            l2 = obj.link2.length;
            J =  [ -l1*obj.S1 - l2*obj.S12, -l2*obj.S12;
                    l1*obj.C1 + l2*obj.C12,  l2*obj.C12 ];
        end
        
        function [Qx, Qy] = Hessian (obj)
            l1 = obj.link1.length;
            l2 = obj.link2.length;
            Qx = [ -l1*obj.C1-l2*obj.C12, -l2*obj.C12;
                             -l2*obj.C12, -l2*obj.C12 ];
            Qy = [ -l1*obj.S1-l2*obj.S12, -l2*obj.S12;
                             -l2*obj.S12, -l2*obj.S12 ];
        end
        
        function [ mat, vec ] = inertia_matrix_and_torque_vector (obj)
            l1 = obj.link1.length; lc1 = obj.link1.length_center;
            m1 = obj.link1.mass;    J1 = obj.link1.inertia_of_moment;
            l2 = obj.link2.length; lc2 = obj.link2.length_center;
            m2 = obj.link2.mass;    J2 = obj.link2.inertia_of_moment;
            mat = zeros(2,2);
            mat(1,1) = J1 + m1*lc1^2 + J2 + m2*(l1^2+lc2^2+2*l1*lc2*obj.C2);
            mat(2,2) = J2 + m2*lc2*2;
            mat(1,2) = J2 + m2*(lc2^2 + l1*lc2*obj.C2);
            mat(2,1) = mat(1,2);
            
            om1 = obj.omega1; om2 = obj.omega2;
            h12 = m2*l1*lc2*obj.S2;
            G1  = (m1*lc1 + m2*l1)*[ -obj.S1; obj.C1 ]' * obj.gravity;
            G2 =          m2*lc2*[ -obj.S12; obj.C12 ]' * obj.gravity;
            vec = [  h12*om2^2 + 2*h12*om1*om2 + G1 + G2; ...
                    -h12*om1^2 + G2 ];
        end
    end
end