% Dynamic deformation of an elastic square object (10&times;10)
% �����`�e�����̂̓��I�ȕό` (10&times;10)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
[points, triangles] = rectangular_object(m, n, width, height);

% E = 1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[lambda, mu] = Lame_constants(Young, nu);
[lambda_vis, mu_vis] = Lame_constants(c, nu);

npoints = size(points,2);
ntriangles = size(triangles,1);
elastic = Body(npoints, points, ntriangles, triangles, thickness);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

tp = 0.5; vpush = 0.8*(height/3)/tp;
th = 0.5;
tf = 2.0;

alpha = 1e+6;

index_floor = 1:10;
index_push = 94:97;

% pushing top region
% �㕔�������Ă���
A = elastic.constraint_matrix([index_floor, index_push]);
b0 = zeros(2*(10+4),1);
b1 = [ zeros(2*10,1); 0; -vpush; 0; -vpush; 0; -vpush; 0; -vpush];
interval = [0, tp];
qinit = zeros(4*npoints,1);
square_object_push = @(t,q) square_object_constraint_param(t,q, elastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(square_object_push, interval, qinit);

% holding top region
% �㕔��ێ����Ă���
b0 = [ zeros(2*10,1); 0; -vpush*tp; 0; -vpush*tp; 0; -vpush*tp; 0; -vpush*tp ];
b1 = zeros(2*(10+4),1);
interval = [tp, tp+th];
qinit = q_push(end,:);
square_object_hold = @(t,q) square_object_constraint_param(t,q, elastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(square_object_hold, interval, qinit);

% releasing top region
% �㕔�����
A = elastic.constraint_matrix([index_floor]);
b0 = zeros(2*10,1);
b1 = zeros(2*10,1);
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
square_object_free = @(t,q) square_object_constraint_param(t,q, elastic, A,b0,b1, alpha);
[time_free, q_free] = ode15s(square_object_free, interval, qinit);

time = [time_push; time_hold; time_free];
q = [q_push; q_hold; q_free];

figure('position', [0, 0, 400, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

clf;
for t = 0:0.1:tp+th+tf
    fprintf("time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw(disps);
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40]);
    yticks([-10:10:40]);
    pbaspect([1 1 1]);
    grid on;
    filename = strcat('10_10/deform_', num2str(floor(1000*t),'%04d'), '.png');
    saveas(gcf, filename, 'png');
end

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = 0:0.01:tp+th+tf
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw(disps);
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40]);
    yticks([-10:10:40]);
    pbaspect([1 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
    disp(t);
end
M(fr) = getframe(gcf);

v = VideoWriter('square_object_10_10', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

