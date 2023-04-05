body = RigidBody_Cuboid(1, 4, 4, 8);
alpha = 1000; % positive large constant for CSM

ext = @(t) external_torque(t);
rotation_quaternion_ODE = @(t,s) rotation_quaternion_ODE_params (t,s, body, alpha, ext);
%interval = 0:0.0001:20;
tf = 20;
interval = [0,tf];
sinit = [1;0;0;0; 0;0;0;0];
[time, s] = ode45(rotation_quaternion_ODE, interval, sinit);

%making_video (body, time, s, 0.1);
making_tiled_video (body, time, s, 0.1);

function dots = rotation_quaternion_ODE_params (t,s, body, alpha, ext)
    q = s(1:4); dotq = s(5:8);
    ddotq = body.calculate_ddotq (q, dotq, alpha, ext(t));
    dots = [dotq;ddotq];
end

function tau = external_torque(t)
    if t <= 5
        tau = [12.00; 0.00; 0.00];
    elseif t <= 10
        tau = [ 0.00; -12.00; 0.00];
    else
        tau = [0;0;0];
    end
end

function making_video (body, time, s, videostep)

    figure('position', [0, 0, 1200, 1200]);
    set(0,'defaultAxesFontSize',12);
    set(0,'defaultTextFontSize',12);
    
    clf('reset');
    fr = 1;
    clear M;
    ts = time(1);
    tf = time(end);
    for t=ts:videostep:tf
        fprintf("time %f\n", t);
        index = nearest_index(time, t);
        q = s(index, 1:4);
        clf;
        body.draw(zeros(3,1), q);
        xlim([-10,10]); ylim([-10,10]); zlim([-10,10]);
        xlabel('x'); ylabel('y'); zlabel('z');
        pbaspect([1 1 1]); grid on; view([-75, 30]);
        drawnow;
        M(fr) = getframe(gcf);
        fr = fr + 1;
    end
    M(fr) = getframe(gcf);
    
    v = VideoWriter('rotation_quaternion', 'MPEG-4');
    open(v);
    writeVideo(v, M);
    close(v);
end

function making_tiled_video (body, time, s, videostep)

    figure('position', [0, 0, 1600, 1600]);
    set(0,'defaultAxesFontSize',12);
    set(0,'defaultTextFontSize',12);
    
    clf('reset');
    fr = 1;
    clear M;
    ts = time(1);
    tf = time(end);
    for t=ts:videostep:tf
        fprintf("time %f\n", t);
        index = nearest_index(time, t);
        q = s(index, 1:4);
        clf;
        tiledlayout(2,2);
        %
        nexttile;
        body.draw(zeros(3,1), q);
        xlim([-10,10]); ylim([-10,10]); zlim([-10,10]);
        xlabel('x'); ylabel('y'); zlabel('z');
        pbaspect([1 1 1]); grid on; view([-75, 30]);
        %
        nexttile;
        body.draw(zeros(3,1), q);
        xlim([-10,10]); ylim([-10,10]); zlim([-10,10]);
        xlabel('x'); ylabel('y'); zlabel('z');
        pbaspect([1 1 1]); grid on; view([0, 0]);   % x-z
        %
        nexttile;
        body.draw(zeros(3,1), q);
        xlim([-10,10]); ylim([-10,10]); zlim([-10,10]);
        xlabel('x'); ylabel('y'); zlabel('z');
        pbaspect([1 1 1]); grid on; view([-90, 0]);   % y-z
        %
        nexttile;
        body.draw(zeros(3,1), q);
        xlim([-10,10]); ylim([-10,10]); zlim([-10,10]);
        xlabel('x'); ylabel('y'); zlabel('z');
        pbaspect([1 1 1]); grid on; view([0, 90]);   % x-y
        
        drawnow;
        M(fr) = getframe(gcf);
        fr = fr + 1;
    end
    M(fr) = getframe(gcf);
    
    v = VideoWriter('rotation_quaternion_tiled', 'MPEG-4');
    open(v);
    writeVideo(v, M);
    close(v);
end
