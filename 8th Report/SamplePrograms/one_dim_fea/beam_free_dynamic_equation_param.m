function dotq = beam_supported_dynamic_equation_param (t,q, n, M,B,K, force)
% dynamic equation of beam (free motion)
%   M:inertia matrix  B:damping matrix  K:stiffness matrix
%   force:function calculating external force applied to a specifed nodal point

% �r�[���̉^���������i���R�^���j
%   M:�����s��  B:�S���s��  K:�����s��
%   force:�ߓ_�ɍ�p����O�͂��v�Z����֐�

    un = q(1:n);
    vn = q(n+1:2*n);
    
    mat = M;
    vec = -K*un-B*vn + force(t,n,un,vn);
    sol = mat \ vec;
    dotvn = sol(1:n);
    
    dotq = [ vn; dotvn ];
end
