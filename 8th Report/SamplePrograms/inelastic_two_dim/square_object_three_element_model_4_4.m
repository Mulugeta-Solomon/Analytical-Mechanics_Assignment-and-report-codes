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
