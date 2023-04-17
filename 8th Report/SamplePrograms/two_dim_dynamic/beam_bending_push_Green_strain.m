% g, cm, sec

addpath('../two_dim_fea');

width = 10; height = 2; thickness = 1;
m = 11; n = 3;
[ points, triangles ] = rectangular_object( m, n, width, height );

% E = 0.1 MPa; c = 0.004 kPa s; rho = 0.020 g/cm^3
Young = 1.0*1e+6; c = 0.04*1e+3; nu = 0.48; density = 0.020;

[ lambda, mu ] = Lame_constants( Young, nu );
[lambda_vis, mu_vis] = Lame_constants(c, nu);
