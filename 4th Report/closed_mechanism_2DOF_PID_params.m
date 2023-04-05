function dotq = closed_mechanism_2DOF_PID_params (t, q, robot, thetad, Kp, Ki, Kd)
% equation of motion of 2DOF closed mechanism under PID control
    disp(num2str(t,"%8.6f"));
    persistent alpha
    if isempty(alpha)
        alpha = 1000;
    end
    
    theta = q(1:4);
    omega = q(5:8);
    etha = q(9:10);
    
    dottheta = omega;
    
    theta12 = theta(1:2); omega12 = omega(1:2);
    theta34 = theta(3:4); omega34 = omega(3:4);
    theta13 = theta([1,3],:); 
    robot.left_arm = robot.left_arm.joint_angles (theta12, omega12);
    robot.right_arm = robot.right_arm.joint_angles (theta34, omega34);
    
    J12 = robot.left_arm.Jacobian;
    J34 = robot.right_arm.Jacobian;    
    [ Q12x, Q12y ] = robot.left_arm.Hessian;
    [ Q34x, Q34y ] = robot.right_arm.Hessian;
    
    [inertia_matrix_12, torque_vector_12] = robot.left_arm.inertia_matrix_and_torque_vector;
    [inertia_matrix_34, torque_vector_34] = robot.right_arm.inertia_matrix_and_torque_vector;
    
    R = robot.left_arm.tip_point - robot.right_arm.tip_point;
    dotR = J12*omega12 - J34*omega34;

    C = [ omega12'*Q12x*omega12 - omega34'*Q34x*omega34; ...
          omega12'*Q12y*omega12 - omega34'*Q34y*omega34 ];
    C = C + 2*alpha*dotR + (alpha^2)*R;

    dotetha = theta13 - thetad;
    %[theta(1) - thetad(1); theta(3) - thetad(2)];
    %theta - thetad;
    tau_left  = [ -Kp(1)*(theta(1) - thetad(1)) - Ki(1)*etha(1) - Kd(1)*omega(1); 0 ];
    tau_right = [ -Kp(2)*(theta(3) - thetad(2)) - Ki(2)*etha(2) - Kd(2)*omega(3); 0 ];
    
    A = [ inertia_matrix_12,        zeros(2,2), -J12'; ...
                 zeros(2,2), inertia_matrix_34,  J34'; ...
                       -J12,               J34, zeros(2,2) ];
    b = [ torque_vector_12 + tau_left; ...
          torque_vector_34 + tau_right; ...
          C ];
    
    p = A \ b;
    dotomega = p(1:4);
    
    dotq = [dottheta; dotomega; dotetha];
end
