% 平衡式と制約式をまとめた等式制約とする
% 外力 fn は圧力による力であり，定ベクトルではないので
% 等式制約を nonlcon に記述する
% g, cm, sec

dt0 = datetime;
addpath('../three_dim_fea');

% E = 0.1 MPa; c = 0.04 kPa s; rho = 1 g/cm^2
Young = 1.0*1e+6; c = 0.4*1e+3; nu = 0.48; density = 1.00;
[ lambda, mu ] = Lame_constants( Young, nu );

l = 16; m = 16; n = 2; membrane_size = 10; thickness = 0.5;
% p = 2e+4 = 0.02e+6 --> 0.002 MPa = 2 KPa
pressure = 40e+4;
%pressure = 80e+4;