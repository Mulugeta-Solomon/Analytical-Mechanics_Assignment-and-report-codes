function dotq = beam_supported_dynamic_equation_param (t, q, n, M,B,K, constraints, alpha, force)
% dynamic equation of beam (under constraints)
%   M:inertia matrix  B:damping matrix  K:stiffness matrix
%   constraints:constraint matrix  alpha:positive constant for CSM
%   force:function calculating external force applied to a specifed nodal point

% ビームの運動方程式（制約がある場合）
%   M:慣性行列  B:粘性行列  K:剛性行列
%   constraints:制約を表す行列  alpha:制約安定化法の定数
%   force:節点に作用する外力を計算する関数

    un = q(1:n);
    vn = q(n+1:2*n);
    
    dotun = vn;
    
    an = size(constraints,2);
    mat = [ M,  -constraints; ...
           -constraints',  zeros(an,an) ];
    vec = [ -K*un-B*vn + force(t,n,un,vn); ...
            constraints'*(2*alpha*vn + (alpha^2)*un) ];
    sol = mat \ vec;
    dotvn = sol(1:n);
    
    dotq = [ dotun; dotvn ];
end
