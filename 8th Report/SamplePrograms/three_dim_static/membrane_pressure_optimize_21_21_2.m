% 平衡式と制約式をまとめた等式制約とする
% 外力 fn は圧力による力であり，定ベクトルではないので
% 等式制約を nonlcon に記述する
% g, cm, sec

dt0 = datetime;
addpath('../three_dim_fea');

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^2
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[ lambda, mu ] = Lame_constants( Young, nu );

l = 16; m = 16; n = 2; membrane_size = 10; thickness = 0.5;
% p = 2e+4 = 0.02e+6 --> 0.002 MPa = 2 KPa
pressure = 40e+4;
%pressure = 80e+4;

[ points, tetrahedra ] = cuboidal_object( l, m, n, membrane_size, membrane_size, thickness );
npoints = size(points,2);
ntetrahedra = size(tetrahedra,1);
elastic = Body(npoints, points, ntetrahedra, tetrahedra);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.calculate_stiffness_matrix;

disps_natural = zeros(3,npoints);
un_natural = reshape(disps_natural, [3*npoints,1]);

%figure('position', [0, 0, 600, 510]);
figure('position', [0, 0, 400, 340]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

clf;
draw_body(elastic,disps_natural);
saveas(gcf, 'membrane_pressure_optimize_21_21_2_natural.png');

K = elastic.Stiffness_Matrix;
index = sort([ 1:l, (1:m-2)*l+1, (1:m-2)*l+l, (m-1)*l+(1:l) ]);
A = elastic.constraint_matrix(index);
numConstraints = 3*size(index,2);
mat = [ K, -A; -A', zeros(numConstraints,numConstraints) ];

positions = elastic.positional_vectors;
boundary_faces = elastic.BoundaryFaces;
bottom_faces = [];
for p=1:size(boundary_faces,1)
    i=boundary_faces(p,1); j=boundary_faces(p,2); k=boundary_faces(p,3);
    if positions(3,i) <= 0 && positions(3,j) <= 0 && positions(3,k) <= 0
        bottom_faces = [ bottom_faces; boundary_faces(p,:) ];
    end
end