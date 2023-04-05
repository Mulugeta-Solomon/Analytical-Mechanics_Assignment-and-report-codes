% PID control of 2DOF closed mechanism

len = 2.00; radius = 0.05; density = 1;
link1 = Link_Cylinder ( len, radius, density );
link2 = Link_Cylinder ( len, radius, density );
base_left = [-1; 0];
link3 = Link_Cylinder ( len, radius, density );
link4 = Link_Cylinder ( len, radius, density );
base_right = [1; 0];
grav = [0; -9.8];
left_arm  = Open_Mechanism_Two_DOF (link1, link2, base_left,  grav);
right_arm = Open_Mechanism_Two_DOF (link3, link4, base_right, grav);
robot = Closed_Mechanism_Two_DOF (left_arm, right_arm);

thetad = [ pi/3; pi/6 ]; Kp = [ 10; 15 ]; Kd = [ 2.5; 2.5]; Ki = [ 1.5; 5.5];
closed_mechanism_2DOF_PD_ode = @(t,q) closed_mechanism_2DOF_PID_params (t, q, robot, thetad, Kp, Ki, Kd);

interval = [ 0, 10 ];
thetainit = [ pi/2; -pi/6; pi/2; pi/6 ]; omegainit = [ 0; 0; 0; 0 ]; ethainit = [ 0; 0];
qinit = [thetainit; omegainit; ethainit];
[time, q] = ode45(closed_mechanism_2DOF_PD_ode, interval, qinit);

thetad1 = thetad(1)*ones(size(time));
thetad3 = thetad(2)*ones(size(time));
plot(time, q(:,1), 'r-', time, q(:,3), 'b-', ...
     time, thetad1, 'r--', time, thetad3, 'b--');
title('PID angles');
xlabel('time');
ylabel('theta1 and theta3');
legend('theta1','theta3');
grid on;
saveas(gcf, 'closed_mechanism_2DOF_PID_angles.png');
fprintf("pause: press any key\n");
pause;

clf;
sz = size(time);
for k=1:100:sz(1)
    theta = q(k,1:4);
    omega = q(k,5:8);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-4,4]);
    ylim([ 0,4]);
    pbaspect([2 1 1]);
    grid on;
    pause(0.01);
end
fprintf("pause: press any key\n");
pause;

clf;
hold on;
for k=1:100:sz(1)
    theta = q(k,1:4);
    omega = q(k,5:8);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-4,4]);
    ylim([ 0,4]);
    pbaspect([2 1 1]);
    grid on;
    pause(0.01);
end
hold off;
saveas(gcf, 'closed_mechanism_2DOF_PID_draw.png');

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = ts:0.01:te
    clf;
    k = nearest_index(time, t);
    theta = q(k,1:4);
    omega = q(k,5:8);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-4,4]);
    ylim([ 0,4]);
    pbaspect([2 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    M(fr) = getframe(gcf);
    fr = fr + 1;
end
M(fr) = getframe(gcf);

v = VideoWriter('closed_mechanism_2DOF_PID_draw', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);
