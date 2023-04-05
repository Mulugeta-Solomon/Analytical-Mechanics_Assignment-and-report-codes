%2DOF open mechanism driven by external torques

len = 2.00; radius = 0.05; density = 1;
link1 = Link_Cylinder ( len, radius, density );
link2 = Link_Cylinder ( len, radius, density );
base = [0; 0];
grav = [0; -9.8];
robot = Open_Mechanism_Two_DOF (link1, link2, base, grav);

%external_torques = @(t) [ 0; 0 ];
period = 2; omega = 2*pi/period;
external_torques = @(t) [ 0.1*sin(omega*t)+0.2; -0.05*sin(omega*t)-0.01 ];

open_mechanism_2DOF_external_torques_ode = @(t,q) open_mechanism_2DOF_external_torques_params (t, q, robot, external_torques);

interval = [ 0, 10 ];
thetainit = [ 0; 0 ]; omegainit = [ 0; 0 ];
qinit = [thetainit; omegainit];
[time, q] = ode45(open_mechanism_2DOF_external_torques_ode, interval, qinit);

plot(time, q(:,1), 'r-', time, q(:,2), 'b-');
grid on;
saveas(gcf, 'open_mechanism_2DOF_external_torques_angles.png');
fprintf("pause: press any key\n");
pause;

clf;
sz = size(time);
for k=1:1:sz(1)
    theta = q(k,1:2);
    omega = q(k,3:4);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-5,5]);
    ylim([-5,5]);
    pbaspect([1 1 1]);
    grid on;
    pause(0.01);
end
fprintf("pause: press any key\n");
pause;

clf;
hold on;
for k=1:1:sz(1)
    theta = q(k,1:2);
    omega = q(k,3:4);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-5,5]);
    ylim([-5,5]);
    pbaspect([1 1 1]);
    grid on;
    pause(0.01);
end
hold off;
saveas(gcf, 'open_mechanism_2DOF_external_torques_draw.png');

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = ts:0.01:te
    clf;
    k = nearest_index(time, t);
    theta = q(k,1:2);
    omega = q(k,3:4);
    robot = robot.joint_angles (theta, omega);
    robot.draw;
    xlim([-5,5]);
    ylim([-5,5]);
    pbaspect([1 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    M(fr) = getframe(gcf);
    fr = fr + 1;
end
M(fr) = getframe(gcf);

v = VideoWriter('open_mechanism_2DOF_external_torques_draw', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);