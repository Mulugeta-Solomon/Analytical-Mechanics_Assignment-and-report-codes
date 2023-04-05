function dotq = damped_pendulum_Cartesian (t,q)
% equation of motion of simple pendulum (Cartesian)
    global mass; global length; global grav; global alpha; global viscous;
    x = q(1); y = q(2); vx = q(3); vy = q(4);
    
    dotx = vx; doty = vy;
    R = sqrt(x^2+(y-length)^2) - length;
    P = 1/sqrt(x^2+(y-length)^2); 
    Rx = x*P; Ry = (y-length)*P;
    Rxx = P - x^2*P^3; 
    Ryy = P - (y-length)^2*P^3;
    Rxy = -x*(y-length)*P^3;
    C = [vx,vy]*[Rxx, Rxy; Rxy, Ryy]*[vx;vy] ...
        + 2*alpha*[Rx, Ry]*[vx;vy]...
        + alpha^2*R;
    fvx = - viscous*vx;
    fvy = - viscous*vy;
    A = [mass, 0, -Rx; 0, mass, -Ry; -Rx, -Ry, 0];
    b = [ fvx; -mass*grav - fvy; C ];
    s = A \ b;
    dotvx = s(1); dotvy = s(2);
    dotq = [dotx; doty; dotvx; dotvy];
end