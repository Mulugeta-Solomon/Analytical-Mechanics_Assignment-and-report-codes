% Dynamic deformation of an elasto-plastic square object (10&times;10)
% �����`�e�Y�����̂̓��I�ȕό` (10&times;10)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
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

index_floor = 1:10;
index_push = 94:97;

% pushing top region
% �㕔�������Ă���
A = elastoplastic.constraint_matrix([index_floor, index_push]);
b0 = zeros(2*(10+4),1);
b1 = [ zeros(2*10,1); 0; -vpush; 0; -vpush; 0; -vpush; 0; -vpush ];
interval = [0, tp];
qinit = zeros(8*npoints,1);
square_object_push = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(square_object_push, interval, qinit);

% holding top region
% �㕔��ێ����Ă���
b0 = [ zeros(2*10,1); 0; -vpush*tp; 0; -vpush*tp; 0; -vpush*tp; 0; -vpush*tp ];
b1 = zeros(2*(10+4),1);
interval = [tp, tp+th];
qinit = q_push(end,:);
square_object_hold = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(square_object_hold, interval, qinit);


% releasing top region
% �㕔�����
A = elastoplastic.constraint_matrix([index_floor]);
b0 = zeros(2*10,1);
b1 = zeros(2*10,1);
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
square_object_free = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_free, q_free] = ode15s(square_object_free, interval, qinit);
