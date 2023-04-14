function f = floor_force_param ( t, n, un, vn, k, A, Evfloor )
% floor reaction force is applied to point P_k
% 点 P_k に床反力が作用する
    f = zeros(n,1);
    if (un(k) < 0.00)
        f(k) = -Evfloor*un(k)*A;
    end
end
