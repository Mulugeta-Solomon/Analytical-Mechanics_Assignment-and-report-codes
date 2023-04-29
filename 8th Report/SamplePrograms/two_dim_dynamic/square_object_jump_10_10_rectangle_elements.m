% Jumping of an elastic square object (10&times;10) rectangular elements
% �����`�e�����̂̒��� (10&times;10) �����`�v�f
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
[ points, triangles, rectangles ] = rectangular_object(m, n, width, height);
