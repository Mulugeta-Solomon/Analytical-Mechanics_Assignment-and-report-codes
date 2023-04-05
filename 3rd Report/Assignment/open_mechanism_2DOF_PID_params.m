function dotq = open_mechanism_2DOF_PID_params (t, q, robot, thetad, Kp, Ki, Kd)
% equation of motion of 2DOF open mechanism under PD control
    disp(num2str(t,"%8.6f"));
    
    theta = q(1:2);
    omega = q(3:4);
    etha = q(5:6);
    
    robot = robot.joint_angles(theta, omega);
    [inertia_matrix, torque_vector] = robot.inertia_matrix_and_torque_vector;
    
    dottheta = omega;
    dotetha = theta - thetad;
    tau = -Kp.*(theta - thetad) - Ki.*etha - Kd.*omega;
    dotomega = inertia_matrix \ (torque_vector + tau);
    
    dotq = [dottheta; dotomega; dotetha];
end