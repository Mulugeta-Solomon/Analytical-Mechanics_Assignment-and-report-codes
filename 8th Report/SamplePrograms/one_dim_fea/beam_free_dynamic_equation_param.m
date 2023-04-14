function dotq = beam_supported_dynamic_equation_param (t,q, n, M,B,K, force)
% dynamic equation of beam (free motion)
%   M:inertia matrix  B:damping matrix  K:stiffness matrix
%   force:function calculating external force applied to a specifed nodal point

% ビームの運動方程式（自由運動）
%   M:慣性行列  B:粘性行列  K:剛性行列
%   force:節点に作用する外力を計算する関数

    un = q(1:n);
    vn = q(n+1:2*n);
    
    mat = M;
    vec = -K*un-B*vn + force(t,n,un,vn);
    sol = mat \ vec;
    dotvn = sol(1:n);
    
    dotq = [ vn; dotvn ];
end
