% Jumping of an elastic square object (4&times;4) rectangular elements
% �����`�e�����̂̒��� (4&times;4) �����`�v�f
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 4; n = 4;
[ points, triangles, rectangles ] = rectangular_object( m, n, width, height );

% E = 1 MPa; c = 0.04 kPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[lambda, mu] = Lame_constants(Young, nu);
[lambda_vis, mu_vis] = Lame_constants(c, nu);

npoints = size(points,2);
nrectangles = size(rectangles,1);
elastic = Body(npoints, points, [], [], thickness);
elastic = elastic.rectangle_elements(nrectangles, rectangles);
elastic = elastic.mechanical_parameters(density, lambda, mu);
elastic = elastic.viscous_parameters(lambda_vis, mu_vis);
elastic = elastic.calculate_stiffness_matrix;
elastic = elastic.calculate_damping_matrix;
elastic = elastic.calculate_inertia_matrix;

tp = 0.5; vpush = 0.8*(height/3)/tp;
th = 0.5;
tf = 0.5;

alpha = 1e+6;

% pushing top region
% �㕔�������Ă���
A = elastic.constraint_matrix([1,2,3,4,14,15]);
b0 = zeros(2*6,1);
b1 = [ zeros(2*4,1); 0; -vpush; 0; -vpush ];
interval = [0, tp];
qinit = zeros(4*npoints,1);
square_object_push = @(t,q) square_object_constraint_param(t,q, elastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(square_object_push, interval, qinit);

% holding top region
% �㕔��ێ����Ă���
b0 = [ zeros(2*4,1); 0; -vpush*tp; 0; -vpush*tp ];
b1 = zeros(2*6,1);
interval = [tp, tp+th];
qinit = q_push(end,:);
square_object_hold = @(t,q) square_object_constraint_param(t,q, elastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(square_object_hold, interval, qinit);

% releasing all constraints
% ���ׂĂ̐�������
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
square_object_free = @(t,q) square_object_free_param(t,q, elastic);
[time_free, q_free] = ode15s(square_object_free, interval, qinit);

time = [time_push; time_hold; time_free];
q = [q_push; q_hold; q_free];

figure('position', [0, 0, 400, 800]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);
floor_color = [0.85 0.85 0.85];

clf;
for t = 0:0.1:tp+th+tf
    fprintf("time %f\n", t);
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw_individual_rectangles(disps);
    fill([40, 40, -10, -10], [-10, 0, 0, -10], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
    hold off;
    xlim([-10,40]);
    ylim([-10,90]);
    xticks([-10:10:40]);
    yticks([-10:10:90]);
    pbaspect([1 2 1]);
    grid on;
    filename = strcat('jump_4_4_rectangle_elements/deform_', num2str(floor(1000*t),'%04d'), '.png');
    saveas(gcf, filename, 'png');
end

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = 0:0.01:tp+th+tf
    index = nearest_index(time, t);
    disps = reshape(q(index,1:npoints*2), [2,npoints]);
    elastic.draw_individual_rectangles(disps);
    fill([40, 40, -10, -10], [-10, 0, 0, -10], floor_color, 'FaceAlpha', 0.2, 'EdgeColor','none');
    hold off;
    xlim([-10,40]);
    ylim([-10,90]);
    xticks([-10:10:40]);
    yticks([-10:10:90]);
    pbaspect([1 2 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
    disp(t);
end
M(fr) = getframe(gcf);

v = VideoWriter('square_object_jump_4_4_rectangle_elements', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = square_object_constraint_param(t,q, body, A,b0,b1, alpha)
    %disp(t);
    
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
    vec = [ -K*un-B*vn; 2*alpha*(A'*vn-b1)+(alpha^2)*(A'*un-(b0+b1*t)) ];
    sol = coef\vec;
    dotvn = sol(1:2*npoints);
    
    dotq = [dotun; dotvn];
end

function dotq = square_object_free_param(t,q, body)
    %disp(t);
    
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
    
    coef = M;
    vec = -K*un-B*vn + contact_force(t, body, un, vn);
    sol = coef\vec;
    dotvn = sol(1:2*npoints);
    
    dotq = [dotun; dotvn];
end

function f = contact_force(t, body, un, vn)

    persistent npoints xn Kcontact Bcontact;
    if isempty(npoints)
        npoints = body.numNodalPoints;
        for k=1:npoints
            xn = [ xn; body.NodalPoints(k).Coordinates ];
        end
        Kcontact = 0.10e+6;
        Bcontact = 0;
    end
    
    rn = xn + un;
    f = zeros(2*npoints,1);
    for k=1:npoints
        if rn(2*k) < 0
            f(2*k) = -Kcontact*rn(2*k) -Bcontact*vn(2*k);
        end
    end
end