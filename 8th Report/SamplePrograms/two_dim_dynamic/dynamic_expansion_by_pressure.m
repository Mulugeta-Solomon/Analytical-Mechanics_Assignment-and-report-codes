% g, cm, sec

addpath('../two_dim_fea');

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[ lambda, mu ] = Lame_constants( Young, nu );
[lambda_vis, mu_vis] = Lame_constants(c, nu);
% p = 2e+4 = 0.02e+6 --> 0.002 MPa = 2 KPa
%pressure = 2e+4; filename = 'expansion_by_pressure_deformed.png';
%pressure = -2e+4; filename = 'expansion_by_pressure_deformed_negative.png';

mrect = 11; nrect = 2;
[ points, triangles ] = rectangular_object( mrect, nrect, 10, 1 );
thickness = 1;

npoints = size(points,2);
ntriangles = size(triangles,1);
elastic = Body(npoints, points, ntriangles, triangles, thickness);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

index_pressure_area = [mrect:-1:1];

A = elastic.constraint_matrix([1,11,12,22]);
b0 = zeros(2*4,1);
b1 = zeros(2*4,1);

alpha = 1e+6;

tf = 3.00;
interval = [0, tf];
qinit = zeros(4*npoints,1);
expansion = @(t,q) object_constraint_param(t,q, elastic, index_pressure_area, A,b0,b1, alpha);
[time, q] = ode45(expansion, interval, qinit);

figure('position', [0, 0, 600, 420]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);
