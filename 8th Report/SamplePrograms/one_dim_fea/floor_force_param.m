function f = floor_force_param ( t, n, un, vn, k, A, Evfloor )
% floor reaction force is applied to point P_k
% �_ P_k �ɏ����͂���p����
    f = zeros(n,1);
    if (un(k) < 0.00)
        f(k) = -Evfloor*un(k)*A;
    end
end
