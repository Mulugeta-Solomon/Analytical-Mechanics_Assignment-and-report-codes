% Dynamic deformation of an elasto-plastic square object (4&times;4)
% �����`�e�Y�����̂̓��I�ȕό` (4&times;4)
% g, cm, sec

addpath('../two_dim_fea');

width = 30; height = 30; thickness = 1;
m = 4; n = 4;
[points, triangles] = rectangular_object(m, n, width, height);
npoints = size(points,2);
ntriangles = size(triangles,1);
elastoplastic = Body_ThreeElementModel(npoints, points, ntriangles, triangles, thickness);

% E = 1 MPa; c1 = 0.04 kPa s; c2 = 2 MPa s; rho = 1 g/cm^3
Young = 10.0*1e+6; c1 = 0.4*1e+3; c2 = 20*1e+6; nu = 0.48; density = 1.00;
[lambda_hard, mu_hard] = Lame_constants(Young, nu);
[lambda_soft, mu_soft] = Lame_constants(0.2*Young, nu);
[lambda_vis_1, mu_vis_1] = Lame_constants(c1, nu);
[lambda_vis_2, mu_vis_2] = Lame_constants(c2, nu);


elastoplastic = elastoplastic.define_subregion([ 1:6, 13:18 ]);
elastoplastic = elastoplastic.subregion_mechanical_parameters(density, lambda_hard, mu_hard, lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2);
elastoplastic = elastoplastic.subregion_color( [0.85 0.85 0.85] );

elastoplastic = elastoplastic.define_subregion([ 7:12 ]);
elastoplastic = elastoplastic.subregion_mechanical_parameters(density, lambda_soft, mu_soft, lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2);
elastoplastic = elastoplastic.subregion_color( [0.95 0.95 0.95] );

elastoplastic = elastoplastic.calculate_stiffness_matrix;
elastoplastic = elastoplastic.calculate_inertia_matrix;

elastoplastic.draw_individual;
drawnow;

nf = elastoplastic.SubRegions(1).numNodalPoints + elastoplastic.SubRegions(2).numNodalPoints;

tp = 1.0; vpush = 0.8*(height/3)/tp;
th = 1.0;
tf = 5.0;

alpha = 1e+6;

% pushing top region
% �㕔�������Ă���
A = elastoplastic.constraint_matrix([1,2,3,4,14,15]);
b0 = zeros(2*6,1);
b1 = [ zeros(2*4,1); 0; -vpush; 0; -vpush ];
interval = [0, tp];
qinit = zeros(4*npoints+4*nf,1);
square_object_push = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_push, q_push] = ode15s(square_object_push, interval, qinit);

% holding top region
% �㕔��ێ����Ă���
b0 = [ zeros(2*4,1); 0; -vpush*tp; 0; -vpush*tp ];
b1 = zeros(2*6,1);
interval = [tp, tp+th];
qinit = q_push(end,:);
square_object_hold = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_hold, q_hold] = ode15s(square_object_hold, interval, qinit);

% releasing top region
% �㕔�����
A = elastoplastic.constraint_matrix([1,2,3,4]);
b0 = zeros(2*4,1);
b1 = zeros(2*4,1);
interval = [tp+th, tp+th+tf];
qinit = q_hold(end,:);
square_object_free = @(t,q) square_object_constraint_param(t,q, elastoplastic, A,b0,b1, alpha);
[time_free, q_free] = ode15s(square_object_free, interval, qinit);

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
    elastoplastic.draw_individual(disps);
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40]);
    yticks([-10:10:40]);
    pbaspect([1 1 1]);
    grid on;
    filename = strcat('4_4_layered/deform_', num2str(floor(1000*t),'%04d'), '.png');
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
    elastoplastic.draw_individual(disps);
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40]);
    yticks([-10:10:40]);
    pbaspect([1 1 1]);
    title(['time ' num2str(t,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
    disp(t);
end
M(fr) = getframe(gcf);

v = VideoWriter('square_object_three_element_model_4_4_layered', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);

function dotq = square_object_constraint_param(t,q, body, A,b0,b1, alpha)
    disp(t);
    
    persistent npoints nsubregions M ps pe flambda_all_s flambda_all_e fmu_all_s fmu_all_e cl1 cl2 cl3 cm1 cm2 cm3;
    if isempty(npoints)
        npoints = body.numNodalPoints;
        nsubregions = body.numSubRegions;
        M = body.Inertia_Matrix;
        [ ps, pe ] = body.flambda_fmu_location_corresponding_to_subregions;
        flambda_all_s = 4*npoints+1;
        flambda_all_e = flambda_all_s -1 + pe(nsubregions);
        fmu_all_s = flambda_all_e + 1;
        fmu_all_e = fmu_all_s -1 + pe(nsubregions);
        cl1 = zeros(nsubregions, 1); cl2 = zeros(nsubregions, 1); cl3 = zeros(nsubregions, 1);
        cm1 = zeros(nsubregions, 1); cm2 = zeros(nsubregions, 1); cm3 = zeros(nsubregions, 1);
        for p=1:nsubregions
            sub = body.SubRegions(p);
            lambda = sub.lambda; mu = sub.mu;
            lambdav1 = sub.lambda_vis_1; muv1 = sub.mu_vis_1;
            lambdav2 = sub.lambda_vis_2; muv2 = sub.mu_vis_2;
            cl1(p) = lambda/(lambdav1+lambdav2);
            cl2(p) = lambda*lambdav2/(lambdav1+lambdav2);
            cl3(p) = lambdav1*lambdav2/(lambdav1+lambdav2);
            cm1(p) = mu/(muv1+muv2);
            cm2(p) = mu*muv2/(muv1+muv2);
            cm3(p) = muv1*muv2/(muv1+muv2);
        end
    end

    un = q(1:2*npoints);
    vn = q(2*npoints+1:4*npoints);
    flambda_all = q(flambda_all_s:flambda_all_e);
    fmu_all     = q(fmu_all_s:fmu_all_e);
    
    dotun = vn;

    f = body.nodal_forces(flambda_all, fmu_all, nsubregions, ps, pe);
    
    coef = [ M, -A; -A', zeros(size(A,2),size(A,2))];
    vec = [ f; 2*alpha*(A'*vn-b1)+(alpha^2)*(A'*un-(b0+b1*t)) ];
    sol = coef\vec;
    dotvn = sol(1:2*npoints);
    
    dotflambda_all = zeros(pe(nsubregions),1);
    dotfmu_all     = zeros(pe(nsubregions),1);
    
    for p=1:nsubregions
        vp = vn(body.SubRegions(p).suffixes_for_nodal_points);
        dotvp = dotvn(body.SubRegions(p).suffixes_for_nodal_points);
        dotflambda_all(ps(p):pe(p)) = -cl1(p)*flambda_all(ps(p):pe(p)) + cl2(p)*vp + cl3(p)*dotvp;
        dotfmu_all    (ps(p):pe(p)) = -cm1(p)*fmu_all    (ps(p):pe(p)) + cm2(p)*vp + cm3(p)*dotvp;
    end
    
    dotq = [dotun; dotvn; dotflambda_all; dotfmu_all];
end