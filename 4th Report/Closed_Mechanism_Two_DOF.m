classdef Closed_Mechanism_Two_DOF
    properties
        left_arm;
        right_arm;
    end
    methods
        function obj = Closed_Mechanism_Two_DOF (left, right)
            obj.left_arm = left;
            obj.right_arm = right;
        end
        
        function obj = joint_angles (obj, theta, omega)
            obj.left_arm  = obj.left_arm.joint_angles (theta(1:2), omega(1:2));
            obj.right_arm = obj.right_arm.joint_angles(theta(3:4), omega(3:4));
        end
        
        function draw (obj)
            hold on;
            obj.left_arm.draw;
            obj.right_arm.draw;
            hold off;
        end
    end
end
