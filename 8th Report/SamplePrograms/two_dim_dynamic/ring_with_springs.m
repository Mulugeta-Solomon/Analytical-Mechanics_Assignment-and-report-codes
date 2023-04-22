% g, cm, sec

addpath('../two_dim_fea');

period01 = 0.2; period12 = 1.8; period23 = 0.1; period34 = 4.8; period45 = 0.1; period56 = 5.0;
% time1=0.2; time2 = 2.0; time3 = 2.1; time4 = 6.6; time5 = 6.7; time6 = 12.0;
time1 =  0 + period01; time2 = time1 + period12; time3 = time2 + period23; time4 = time3 + period34; time5 = time4 + period45; time6 = time5 + period56;
tinterval = 0.1;
vinterval = 0.01;
[ 0, time1, time2, time3, time4, time5, time6 ]

m = 32; n = 3; router = 5; rinner = 4; thickness = 1;
[points, triangles] = ring_object(m, n, router, rinner);
points = points + [ 0; router ];

npoints = size(points, 2);
ntriangles = size(triangles, 1);
ring = Body(npoints, points, ntriangles, triangles, thickness);

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[ lambda, mu ] = Lame_constants( Young, nu );
[ lambda_vis, mu_vis ] = Lame_constants( c, nu );
ring = ring.mechanical_parameters(density, lambda, mu);
ring = ring.viscous_parameters(lambda_vis, mu_vis);

ring = ring.calculate_stiffness_matrix;
ring = ring.calculate_damping_matrix;
ring = ring.calculate_inertia_matrix;
ring = ring.calculate_coefficient_matrices_for_Green_strain;

% g = 9.8 m/s^2 = 980 cm/s^2
grav = [0; -980];
ring = ring.calculate_gravitational_vector(grav);

mass = 1;

for i=1:8
    springs(i) = Spring(2e+4, 0.0, rinner);
end
connect = [0:7]*(m/8)*n + 1;
%extensional_forces = [ 0; -1e+5; 0; 0; 0; 0; 0; 0 ];
% N = Kg m/s^2 = 10^3 g 10^2 cm/s^2 = 10^5 g cm/s^2

uninit = zeros(2*npoints,1);
vninit = zeros(2*npoints,1);
xcinit = [ 0; router ];
vcinit = [ 0; 0 ];
qinit = [ uninit; vninit; xcinit; vcinit ];


