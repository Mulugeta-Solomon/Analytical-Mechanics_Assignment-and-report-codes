% g, cm, sec

addpath('../two_dim_fea');

width = 10; height = 2; thickness = 1;
m = 11; n = 3;
[ points, triangles ] = rectangular_object( m, n, width, height );

% E = 0.1 MPa; c = 0.004 kPa s; rho = 0.020 g/cm^3
Young = 1.0*1e+6; c = 0.04*1e+3; nu = 0.48; density = 0.020;

[ lambda, mu ] = Lame_constants( Young, nu );
[lambda_vis, mu_vis] = Lame_constants(c, nu);

npoints = size(points,2);
ntriangles = size(triangles,1);
elastic = Body(npoints, points, ntriangles, triangles, thickness);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

tp = 0.3; upush = 0; vpush = -2.5*height/tp;
th = 0.2;
tf = 1.5;

alpha = 1e+6;

% pushing
A = elastic.constraint_matrix([1,12,23, 22]);
b0 = zeros(2*4,1);
b1 = [ zeros(2*3,1); upush; vpush ];
interval = [0,tp];
qinit = zeros(4*npoints,1);
f_beam_bending_push = @(t,q) beam_bending_push_Green_param(t,q, elastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(f_beam_bending_push, interval, qinit);

% holding
b0 = [ zeros(2*3,1); upush*tp; vpush*tp ];
b1 = zeros(2*4,1);
interval = [tp, tp+th];
qinit = q_push(end,:);
f_beam_bending_push = @(t,q) beam_bending_push_Green_param(t,q, elastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(f_beam_bending_push, interval, qinit);

% free
A = elastic.constraint_matrix([1,12,23]);
b0 = zeros(2*3,1);
b1 = zeros(2*3,1);
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
f_beam_bending_push = @(t,q) beam_bending_push_Green_param(t,q, elastic, A,b0,b1, alpha);
[time_free, q_free] = ode15s(f_beam_bending_push, interval, qinit);

time = [time_push; time_hold; time_free];
q = [q_push; q_hold; q_free];

figure('position', [0, 0, 400, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

clf;
for t = 0:0.1:tp+th+tf
    fprintf("time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw(disps);
    hold off;
    xlim([0,12]);
    ylim([-6,6]);
    pbaspect([1 1 1]);
    grid on;
    filename = strcat('beam_bending_push_Green_strain/deform_', num2str(floor(1000*t),'%04d'), '.png');
    saveas(gcf, filename, 'png');
end

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = 0:0.01:tp+th+tf
    fprintf("video time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw(disps);
    hold off;
    xlim([0,12]);
    ylim([-6,6]);
    pbaspect([1 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
end
M(fr) = getframe(gcf);

v = VideoWriter('beam_bending_push_Green_strain', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = beam_bending_push_Green_param(t,q, body, A,b0,b1, alpha)
    %disp(t);
    disp(num2str(t,"%6.4f"));
    
    persistent npoints M B K;
    if isempty(npoints)
        npoints = body.numNodalPoints;
        M = body.Inertia_Matrix;
        B = body.Damping_Matrix;
        K = body.Stiffness_Matrix;
    end
    
    un = q(1:2*npoints);
    vn = q(2*npoints+1:4*npoints);
    
    dotun = vn;
    
    coef = [ M, -A; -A', zeros(size(A,2),size(A,2))];
    forces = body.nodal_forces_Green_strain(reshape(un, [2,npoints]));
    vec = [ forces-B*vn; 2*alpha*(A'*vn-b1)+(alpha^2)*(A'*un-(b0+b1*t)) ];
    sol = coef\vec;
    dotvn = sol(1:2*npoints);
    
    dotq = [dotun; dotvn];
end
