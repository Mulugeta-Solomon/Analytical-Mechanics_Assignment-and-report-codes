% g, cm, sec

dt0 = datetime;
addpath('../two_dim_fea');

m = 16; n = 3; router = 4; rinner = 2; thickness = 1;
[points, triangles] = ring_object(m, n, router, rinner);
points = points + [ 0; router ];

npoints = size(points, 2);
ntriangles = size(triangles, 1);
ring = Body(npoints, points, ntriangles, triangles, thickness);
