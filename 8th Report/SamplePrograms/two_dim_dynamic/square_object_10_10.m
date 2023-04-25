% Dynamic deformation of an elastic square object (10&times;10)
% �����`�e�����̂̓��I�ȕό` (10&times;10)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 10; n = 10;
[points, triangles] = rectangular_object(m, n, width, height);
