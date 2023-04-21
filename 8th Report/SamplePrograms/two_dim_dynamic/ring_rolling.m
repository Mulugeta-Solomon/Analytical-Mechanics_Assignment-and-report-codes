% g, cm, sec

dt0 = datetime;
addpath('../two_dim_fea');

m = 16; n = 3; router = 4; rinner = 2; thickness = 1;
[points, triangles] = ring_object(m, n, router, rinner);
points = points + [ 0; router ];

npoints = size(points, 2);
ntriangles = size(triangles, 1);
ring = Body(npoints, points, ntriangles, triangles, thickness);

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
density_heavy = 10*density;
[ lambda, mu ] = Lame_constants( Young, nu );
[ lambda_vis, mu_vis ] = Lame_constants( c, nu );

ring = ring.define_subregion([1:32]);
ring = ring.subregion_mechanical_parameters(density_heavy, lambda, mu);
ring = ring.subregion_viscous_parameters(lambda_vis, mu_vis);
ring = ring.subregion_color( [0.85 0.85 0.85] );

ring = ring.define_subregion([33:64]);
ring = ring.subregion_mechanical_parameters(density, lambda, mu);
ring = ring.subregion_viscous_parameters(lambda_vis, mu_vis);
ring = ring.subregion_color( [0.95 0.95 0.95] );

figure('position', [0, 0, 600, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);
floor_color = [0.85 0.85 0.85];

clf;
ring.draw_individual;
fill([12, 12, -6, -6], [-2, 0, 0, -2], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
xlim([-6,12]); ylim([-2,10]);
pbaspect([1.5 1 1]);
grid on;
drawnow;

ring = ring.calculate_stiffness_matrix;
ring = ring.calculate_damping_matrix;
ring = ring.calculate_inertia_matrix;
ring = ring.calculate_coefficient_matrices_for_Green_strain;

% g = 9.8 m/s^2 = 980 cm/s^2
grav = [0; -980];
ring = ring.calculate_gravitational_vector(grav);

tf = 5.0;
interval = [0, tf];
ring_free = @(t,q) ring_free_param(t,q, ring);

uninit = zeros(2*npoints,1);
vninit = zeros(2*npoints,1);
qinit = [uninit; vninit];
%[time, q] = ode15s(tube_free, interval, qinit);
[time, q] = ode45(ring_free, interval, qinit);

dt1 = datetime;
dt01 = between(dt0, dt1);
dt01

