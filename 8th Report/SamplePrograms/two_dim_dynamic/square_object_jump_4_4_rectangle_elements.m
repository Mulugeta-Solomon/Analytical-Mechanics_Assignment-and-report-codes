% Jumping of an elastic square object (4&times;4) rectangular elements
% �����`�e�����̂̒��� (4&times;4) �����`�v�f
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 4; n = 4;
[ points, triangles, rectangles ] = rectangular_object( m, n, width, height );
