% Calculation based on Green strain
% �O���[���c�݂ɂ��v�Z

% g, cm, sec

addpath('../three_dim_fea');

l = 2; m = 2; n = 5;
cube_size = 1;
rotation_angle = 20*2*pi/360;   % 20 degrees

[ points, tetrahedra ] = cuboidal_object( l, m, n, cube_size*(l-1), cube_size*(m-1), cube_size*(n-1) );

% E = 0.1 MPa; rho = 0.020 g/cm^2
Young = 1.0*1e+6; nu = 0.48; density = 0.020;
% 10.0 N --> fpush = 10.0*1e+5;
%force = [0; -1.2*1e+5];

[ lambda, mu ] = Lame_constants( Young, nu );

npoints = size(points,2);
ntetrahedra = size(tetrahedra,1);
elastic = Body(npoints, points, ntetrahedra, tetrahedra);
elastic = elastic.mechanical_parameters(density, lambda, mu);

disps_natural = zeros(3,npoints);
un_natural = reshape(disps_natural, [3*npoints,1]);