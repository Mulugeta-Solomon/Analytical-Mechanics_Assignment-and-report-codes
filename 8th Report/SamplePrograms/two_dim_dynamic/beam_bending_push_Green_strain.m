% g, cm, sec

addpath('../two_dim_fea');

width = 10; height = 2; thickness = 1;
m = 11; n = 3;
[ points, triangles ] = rectangular_object( m, n, width, height );
