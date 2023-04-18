% g, cm, sec

addpath('../two_dim_fea');

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[ lambda, mu ] = Lame_constants( Young, nu );
[lambda_vis, mu_vis] = Lame_constants(c, nu);
% p = 2e+4 = 0.02e+6 --> 0.002 MPa = 2 KPa
%pressure = 2e+4; filename = 'expansion_by_pressure_deformed.png';
%pressure = -2e+4; filename = 'expansion_by_pressure_deformed_negative.png';
