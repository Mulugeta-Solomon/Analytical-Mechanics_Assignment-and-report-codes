% Jumping of an elastic square object (10&times;10)
% �����`�e�����̂̒��� (10&times;10)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
[points, triangles] = rectangular_object(m, n, width, height);


% E = 1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[lambda, mu] = Lame_constants(Young, nu);
[lambda_vis, mu_vis] = Lame_constants(c, nu);
