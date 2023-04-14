function dotq = beam_supported_dynamic_equation_param (t, q, n, M,B,K, constraints, alpha, force)
% dynamic equation of beam (under constraints)
%   M:inertia matrix  B:damping matrix  K:stiffness matrix
%   constraints:constraint matrix  alpha:positive constant for CSM
%   force:function calculating external force applied to a specifed nodal point

% �r�[���̉^���������i���񂪂���ꍇ�j
%   M:�����s��  B:�S���s��  K:�����s��
%   constraints:�����\���s��  alpha:������艻�@�̒萔
%   force:�ߓ_�ɍ�p����O�͂��v�Z����֐�

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
