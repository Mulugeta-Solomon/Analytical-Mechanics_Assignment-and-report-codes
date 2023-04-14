function f = pushing_force_param ( t, n, un, vn, k, fpush, tpush )
% external force is applied to point P_k
% �_ P_k �ɊO�͂���p����
    f = zeros(n,1);
    
    % time interval [0,tpush] : downward force -fpush
    % ���ԋ�� [0,tpush] : P_k �ɉ������O�� -fpush
    f(k) = -fpush * (t <= tpush);
end
