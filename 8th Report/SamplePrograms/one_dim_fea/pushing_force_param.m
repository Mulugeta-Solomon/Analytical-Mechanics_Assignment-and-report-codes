function f = pushing_force_param ( t, n, un, vn, k, fpush, tpush )
% external force is applied to point P_k
% 点 P_k に外力が作用する
    f = zeros(n,1);
    
    % time interval [0,tpush] : downward force -fpush
    % 時間区間 [0,tpush] : P_k に下向き外力 -fpush
    f(k) = -fpush * (t <= tpush);
end
