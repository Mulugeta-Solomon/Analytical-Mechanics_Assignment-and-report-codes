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

clf;
for t = 0:0.1:tf
    fprintf("time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    ring.draw_individual(disps);
    fill([12, 12, -6, -6], [-2, 0, 0, -2], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
    hold off;
    xlim([-6,12]); ylim([-2,10]); xticks([-6:2:12]);
    pbaspect([1.5 1 1]);
    grid on;
    filename = strcat('ring_rolling/deform_', num2str(floor(1000*t),'%04d'), '.png');
    saveas(gcf, filename, 'png');
end

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = 0:0.01:tf
    fprintf("video time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    ring.draw_individual(disps);
    fill([12, 12, -6, -6], [-2, 0, 0, -2], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
    hold off;
    xlim([-6,12]); ylim([-2,10]); xticks([-6:2:12]);
    pbaspect([1.5 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
end
M(fr) = getframe(gcf);

v = VideoWriter('ring_rolling', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = ring_free_param(t,q, body)
    %disp(t);
    disp(num2str(t,"%8.6f"));
    
    persistent npoints xn thickness M B K gravitational_force;
    persistent Kcontact Bcontact friction_damping;
    if isempty(npoints)
        npoints = body.numNodalPoints;
        for k=1:npoints
            xn = [ xn; body.NodalPoints(k).Coordinates ];
        end
        thickness = body.Thickness;
        M = body.Inertia_Matrix;
        B = body.Damping_Matrix;
        K = body.Stiffness_Matrix;
        gravitational_force = body.Gravitational_Vector;
        Kcontact = 1e+6; % N/m = kg/s^2 = 10^3 g/s^2
        Bcontact = 0;
        friction_damping = 2000; % Ns/m = 10^3 g/s
    end
    
    un = q(1:2*npoints);
    vn = q(2*npoints+1:4*npoints);
    disps = reshape(un, [2,npoints]);
    
    % Cauchy strain
    %elastic_force = -K*un -B*vn;
    % Green strain
    elastic_force = body.nodal_forces_Green_strain(disps) -B*vn;
    
    rn = xn + un;
    contact_force = zeros(2*npoints,1);
    for k=1:npoints
        if rn(2*k) < 0
            contact_force(2*k-1) = -friction_damping*vn(2*k-1);
            contact_force(2*k) = -Kcontact*rn(2*k) -Bcontact*vn(2*k);            
        end
    end
    
    f = elastic_force + gravitational_force + contact_force;
    
    dotun = vn;
    
    coef = M;
    vec = f;
    sol = coef\vec;
    dotvn = sol(1:2*npoints);
    
    dotq = [dotun; dotvn];
end
