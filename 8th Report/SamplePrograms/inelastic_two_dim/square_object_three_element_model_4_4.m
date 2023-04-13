% Dynamic deformation of an elasto-plastic square object (4&times;4)
% �����`�e�Y�����̂̓��I�ȕό` (4&times;4)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 4; n = 4;
[points, triangles] = rectangular_object(m, n, width, height);
npoints = size(points,2);
ntriangles = size(triangles,1);
elastoplastic = Body_ThreeElementModel(npoints, points, ntriangles, triangles, thickness);

% E = 1 MPa; c1 = 0.04 kPa s; c2 = 2 MPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c1 = 0.4*1e+3; c2 = 20*1e+6; nu = 0.48; density = 1.00;
[lambda, mu] = Lame_constants(Young, nu);
[lambda_vis_1, mu_vis_1] = Lame_constants(c1, nu);
[lambda_vis_2, mu_vis_2] = Lame_constants(c2, nu);

elastoplastic = elastoplastic.mechanical_parameters(density, lambda, mu, lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2);
elastoplastic = elastoplastic.calculate_stiffness_matrix;
elastoplastic = elastoplastic.calculate_inertia_matrix;

tp = 1.0; vpush = 0.8*(height/3)/tp;
th = 1.0;
tf = 5.0;

alpha = 1e+6;

% pushing top region
% �㕔�������Ă���
A = elastoplastic.constraint_matrix([1,2,3,4,14,15]);
b0 = zeros(2*6,1);
b1 = [ zeros(2*4,1); 0; -vpush; 0; -vpush ];
interval = [0, tp];
qinit = zeros(8*npoints,1);
square_object_push = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(square_object_push, interval, qinit);

% holding top region
% �㕔��ێ����Ă���
b0 = [ zeros(2*4,1); 0; -vpush*tp; 0; -vpush*tp ];
b1 = zeros(2*6,1);
interval = [tp, tp+th];
qinit = q_push(end,:);
square_object_hold = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(square_object_hold, interval, qinit);


% releasing top region
% �㕔�����
A = elastoplastic.constraint_matrix([1,2,3,4]);
b0 = zeros(2*4,1);
b1 = zeros(2*4,1);
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
square_object_free = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
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
    elastoplastic.draw(disps);
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40]);
    yticks([-10:10:40]);
    pbaspect([1 1 1]);
    grid on;
    filename = strcat('4_4/deform_', num2str(floor(1000*t),'%04d'), '.png');
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
    elastoplastic.draw(disps);
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

v = VideoWriter('square_object_three_element_model_4_4', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = square_object_constraint_param(t,q, body, A,b0,b1, alpha)
    disp(t);
    
    persistent npoints lambda mu lambdav1 muv1 lambdav2 muv2 M Jlambda Jmu;
    if isempty(npoints)
        npoints = body.numNodalPoints;
        lambda = body.lambda; mu = body.mu;
        lambdav1 = body.lambda_vis_1; muv1 = body.mu_vis_1;
        lambdav2 = body.lambda_vis_2; muv2 = body.mu_vis_2;
        M = body.Inertia_Matrix;
        Jlambda = body.J_lambda;
        Jmu     = body.J_mu;
    end