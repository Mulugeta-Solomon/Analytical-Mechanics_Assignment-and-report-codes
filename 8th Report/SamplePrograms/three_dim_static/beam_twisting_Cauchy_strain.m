% Calculation based on Cauchy strain
% �R�[�V�[�c�݂ɂ��v�Z

% g, cm, sec

addpath('../three_dim_fea');

l = 2; m = 2; n = 5;
cube_size = 1;
rotation_angle = 20*2*pi/360;   % 20 degrees


[points, tetrahedra ] = cuboidal_object( l, m, n, cube_size*(l-1), cube_size*(m-1), cube_size*(n-1) );

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

internal_energy = @(un) internal_energy_params( un, npoints, elastic );

% fix point 1, 2, 3, 4, lmn-3, lmn-2, lmn-1, lmn
% a set of constraints: A' un = b (each column specifies one constraint)
A = zeros(3*npoints, 3*8); b = zeros(3*8,1);
k = 1; A(3*k-2, 1) = 1; A(3*k-1, 2) = 1; A(3*k, 3) = 1;
k = 2; A(3*k-2, 4) = 1; A(3*k-1, 5) = 1; A(3*k, 6) = 1;
k = 3; A(3*k-2, 7) = 1; A(3*k-1, 8) = 1; A(3*k, 9) = 1;
k = 4; A(3*k-2,10) = 1; A(3*k-1,11) = 1; A(3*k,12) = 1;
p = l*m*n;
disp_corner = @(r, theta) r*[ cos(theta+rotation_angle)-cos(theta); sin(theta+rotation_angle)-sin(theta); 0 ];
k = p-3; A(3*k-2,13) = 1;  A(3*k-1,14) = 1; A(3*k,15) = 1; b(13:15) = disp_corner(cube_size/sqrt(2), (-3/4)*pi);
k = p-2; A(3*k-2,16) = 1;  A(3*k-1,17) = 1; A(3*k,18) = 1; b(16:18) = disp_corner(cube_size/sqrt(2), (-1/4)*pi);
k = p-1; A(3*k-2,19) = 1;  A(3*k-1,20) = 1; A(3*k,21) = 1; b(19:21) = disp_corner(cube_size/sqrt(2), ( 3/4)*pi);
k = p-0; A(3*k-2,22) = 1;  A(3*k-1,23) = 1; A(3*k,24) = 1; b(22:24) = disp_corner(cube_size/sqrt(2), ( 1/4)*pi);

options = optimset('Display','iter', 'MaxFunEvals',3*10^5);
[un_Cauchy, Emin_Cauchy] = fmincon(internal_energy, un_natural, [], [], A', b, [], [], [], options);

figure('position', [0, 0, 400, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

clf;
draw_body(elastic, reshape(un_natural, [3,npoints]));
saveas(gcf, 'beam_twisting_Cauchy_strain_natural.png');

clf;
draw_body(elastic, reshape(un_Cauchy, [3,npoints]));
saveas(gcf, 'beam_twisting_Cauchy_strain_deformed.png');

function energy = internal_energy_params( un, npoints, elastic )
    disps = reshape(un, [3,npoints]);
    energy = elastic.total_strain_potential_energy(disps);
end

function draw_body (body, disps)
    body.draw(disps);
    xlim([-1,2]); xlabel('x');
    ylim([-1,2]); ylabel('y');
    zlim([-1,5]); zlabel('z'); zticks(-1:5);
    pbaspect([1 1 2]);
    grid on;
    view([-75, 30]);
end