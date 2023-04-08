% Jumping of an elastic square object (10&times;10) rectangular elements
% �����`�e�����̂̒��� (10&times;10) �����`�v�f
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
[ points, triangles, rectangles ] = rectangular_object(m, n, width, height);

% E = 1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[lambda, mu] = Lame_constants(Young, nu);
[lambda_vis, mu_vis] = Lame_constants(c, nu);


npoints = size(points,2);
nrectangles = size(rectangles,1);
elastic = Body(npoints, points, [], [], thickness);
elastic = elastic.rectangle_elements(nrectangles, rectangles);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

tp = 0.5; vpush = 0.8*(height/3)/tp;
th = 0.5;
tf = 0.5;

alpha = 1e+6;

index_floor = 1:10;
index_push = 94:97;

