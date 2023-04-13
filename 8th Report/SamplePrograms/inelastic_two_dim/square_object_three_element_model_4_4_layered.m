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
[lambda_hard, mu_hard] = Lame_constants(Young, nu);
[lambda_soft, mu_soft] = Lame_constants(0.2*Young, nu);
[lambda_vis_1, mu_vis_1] = Lame_constants(c1, nu);
[lambda_vis_2, mu_vis_2] = Lame_constants(c2, nu);


elastoplastic = elastoplastic.define_subregion([ 1:6, 13:18 ]);
elastoplastic = elastoplastic.subregion_mechanical_parameters(density, lambda_hard, mu_hard, lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2);
elastoplastic = elastoplastic.subregion_color( [0.85 0.85 0.85] );
