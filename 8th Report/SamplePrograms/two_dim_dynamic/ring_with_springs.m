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

%load('../two_dim_static/ring_and_springs_deformed_grav.mat', 'un_deformed', 'xc_deformed');
%uninit = un_deformed;
%xcinit = xc_deformed;
%qinit = [ uninit; vninit; xcinit; vcinit ];

figure('position', [0, 0, 800, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);
floor_color = [0.85 0.85 0.85];

clf;
%ring.draw;
disps = reshape(uninit, [2,npoints]);
ring.draw_individual(disps);
draw_mass_and_springs(ring, xcinit, connect, disps);
fill([22, 22, -6, -6], [-2, 0, 0, -2], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
xlim([-6,22]); ylim([-2,12]);
pbaspect([2 1 1]);
grid on;
drawnow;

interval = [0, time1];
extensional_forces = [ 0; 0; 0; 0; 0; 0; 0; 0 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t1, q1] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t1, q1, tinterval, ring, connect, floor_color );
qinit = q1(end,:);
M01 = make_video_clip( gcf, t1, q1, vinterval, ring, connect, floor_color );
save('ring_and_springs_period01.mat', 't1', 'q1', 'M01', '-v7.3');
clear t1 q1;

interval = [time1, time2];
extensional_forces = [ 0; -1e+5; 0; 0; 0; -1e+5; 0; 0 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t2, q2] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t2, q2, tinterval, ring, connect, floor_color );
qinit = q2(end,:);
M12 = make_video_clip( gcf, t2, q2, vinterval, ring, connect, floor_color );
save('ring_and_springs_period12.mat', 't2', 'q2', 'M12', '-v7.3');
clear t2 q2;

interval = [time2, time3];
extensional_forces = [ 0; 0; 0; 0; 0; 0; 0; 0 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t3, q3] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t3, q3, tinterval, ring, connect, floor_color );
qinit = q3(end,:);
M23 = make_video_clip( gcf, t3, q3, vinterval, ring, connect, floor_color );
save('ring_and_springs_period23.mat', 't3', 'q3', 'M23', '-v7.3');
clear t3 q3;

interval = [time3, time4];
extensional_forces = [ 0; 0; -1e+5; 0; 0; 0; -1e+5; 0 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t4, q4] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t4, q4, tinterval, ring, connect, floor_color );
qinit = q4(end,:);
M34 = make_video_clip( gcf, t4, q4, vinterval, ring, connect, floor_color );
save('ring_and_springs_period34.mat', 't4', 'q4', 'M34', '-v7.3');
clear t4 q4;

interval = [time4, time5];
extensional_forces = [ 0; 0; 0; 0; 0; 0; 0; 0 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t5, q5] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t5, q5, tinterval, ring, connect, floor_color );
qinit = q5(end,:);
M45 = make_video_clip( gcf, t5, q5, vinterval, ring, connect, floor_color );
save('ring_and_springs_period45.mat', 't5', 'q5', 'M45', '-v7.3');
clear t5 q5;

interval = [time5, time6];
extensional_forces = [ 0; 0; 0; -1e+5; 0; 0; 0; -1e+5 ];
ring_free = @(t,q) ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces);
[t6, q6] = ode45(ring_free, interval, qinit);
draw_ring_and_springs ( gcf, t6, q6, tinterval, ring, connect, floor_color );
qinit = q6(end,:);
M56 = make_video_clip( gcf, t6, q6, vinterval, ring, connect, floor_color );
save('ring_and_springs_period56.mat', 't6', 'q6', 'M56', '-v7.3');
clear t6 q6;

% video clip
M = [ M01, M12(2:end), M23(2:end), M34(2:end), M45(2:end), M56(2:end) ];
v = VideoWriter('ring_with_springs', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = ring_free_param(t,q, ring, grav, mass, springs, connect, extensional_forces)
    %disp(t);
    disp(num2str(t,"%8.6f"));
    
    persistent npoints xn thickness M B K gravitational_force;
    persistent Kcontact Bcontact friction_damping;
    if isempty(npoints)
        npoints = ring.numNodalPoints;
        for k=1:npoints
            xn = [ xn; ring.NodalPoints(k).Coordinates ];
        end
        thickness = ring.Thickness;
        M = ring.Inertia_Matrix;
        B = ring.Damping_Matrix;
        K = ring.Stiffness_Matrix;
        gravitational_force = ring.Gravitational_Vector;
        Kcontact = 1e+6; % N/m = kg/s^2 = 10^3 g/s^2
        Bcontact = 0;
        friction_damping = 20000; % Ns/m = 10^3 g/s
    end
    un = q(1:2*npoints);
    vn = q(2*npoints+1:4*npoints);
    xc = q(4*npoints+1:4*npoints+2);
    vc = q(4*npoints+3:4*npoints+4);

    disps = reshape(un, [2,npoints]);
    
    % Cauchy strain
    %elastic_force = -K*un -B*vn;
    % Green strain
    elastic_force = ring.nodal_forces_Green_strain(disps) -B*vn;
    
    rn = xn + un;
    contact_force = zeros(2*npoints,1);
    for k=1:npoints
        if rn(2*k) < 0
            contact_force(2*k-1) = -friction_damping*vn(2*k-1);
            contact_force(2*k) = -Kcontact*rn(2*k) -Bcontact*vn(2*k);            
        end
    end


