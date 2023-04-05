% solve the equation of motion of damped pendulum (Cartesian)

global mass; global length; global grav;
mass = 0.01; length = 2.0; grav = 9.8;
global alpha;
alpha = 1000;
global viscous;
viscous = 0.01;

interval = [0, 10];
qinit = [length*sin(pi/3); length*(1-cos(pi/3)); 0; 0];
[time, q] = ode45(@damped_pendulum_Cartesian, interval, qinit);

% time - x and y
plot(time,q(:,1),'-', time,q(:,2),'--');
title('time and position');
xlabel('time'); 
ylabel('position X and Y');
legend('X','Y');
saveas(gcf, 'damped_pendulum_Cartesian_x_y.png');
figure;

% time - vx and vy
plot(time,q(:,3),'-', time,q(:,4),'--');
title('time and velocity');
xlabel('time');
ylabel('velocity along X and Y');
legend('vx','vy');
saveas(gcf, 'damped_pendulum_Cartesian_vx_vy.png');
figure;

% x - y
plot(q(:,1), q(:,2));
title('damped pendulum Cartesian path')
xlabel('position along x'); 
ylabel('position along y');
saveas(gcf, 'damped_pendulum_Cartesian_path.png');
figure;

% computed pendulum angle
angle = atan2(q(:,1), length-q(:,2));
plot(time, angle);
title('damped pendulum Cartesian computed angle');
xlabel('time');
ylabel('angle');
saveas(gcf, 'damped_pendulum_Cartesian_computed_angle.png');
figure;

% constraint
R = sqrt(q(:,1).^2 + (q(:,2)-length).^2)-length;
plot(time, R);
title('damped pendulum cartesian R');
xlabel('time');
ylabel('constraint R');
saveas(gcf, 'damped_pendulum_Cartesian_R.png');
figure;