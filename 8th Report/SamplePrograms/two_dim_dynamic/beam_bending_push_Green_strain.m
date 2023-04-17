% g, cm, sec

addpath('../two_dim_fea');

width = 10; height = 2; thickness = 1;
m = 11; n = 3;
[ points, triangles ] = rectangular_object( m, n, width, height );

% E = 0.1 MPa; c = 0.004 kPa s; rho = 0.020 g/cm^3
Young = 1.0*1e+6; c = 0.04*1e+3; nu = 0.48; density = 0.020;

[ lambda, mu ] = Lame_constants( Young, nu );
[lambda_vis, mu_vis] = Lame_constants(c, nu);

npoints = size(points,2);
ntriangles = size(triangles,1);
elastic = Body(npoints, points, ntriangles, triangles, thickness);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

tp = 0.3; upush = 0; vpush = -2.5*height/tp;
th = 0.2;
tf = 1.5;

alpha = 1e+6;

% pushing
A = elastic.constraint_matrix([1,12,23, 22]);
b0 = zeros(2*4,1);
b1 = [ zeros(2*3,1); upush; vpush ];
interval = [0,tp];
qinit = zeros(4*npoints,1);
f_beam_bending_push = @(t,q) beam_bending_push_Green_param(t,q, elastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(f_beam_bending_push, interval, qinit);

% holding
b0 = [ zeros(2*3,1); upush*tp; vpush*tp ];
b1 = zeros(2*4,1);
interval = [tp, tp+th];
qinit = q_push(end,:);
f_beam_bending_push = @(t,q) beam_bending_push_Green_param(t,q, elastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(f_beam_bending_push, interval, qinit);
